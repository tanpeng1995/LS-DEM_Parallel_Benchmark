#ifndef MPI_COLLISION_RESOLUTION_HPP_
#define MPI_COLLISION_RESOLUTION_HPP_
#include "Grain3d.hpp"
#include "definitions.hpp"
#include "utilities.hpp"
#include <algorithm>

struct ContactInfo;
struct WrappedContact{
public:
  ContactInfo _contact;
  Matrix3d     _masterI;
  Matrix3d     _slaveI;
  Matrix3d     _matA;
  Matrix3d     _matB;
  Matrix3d     _colMat;
  double       _mass;

  WrappedContact(const ContactInfo & contact){
    _contact = contact;
    Matrix3d Identity;
    int masterID = contact._master;
    int slaveID = contact._slave;
    Identity << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    double slaveMass = grainsWorld[slaveID].getMass();
    double masterMass = grainsWorld[masterID].getMass();
    if(massScaling) masterMass = masterMass < slaveMass / 10. ? slaveMass / 10. : masterMass;
    _masterI = grainsWorld[masterID].getGlobalI();
    _slaveI  = grainsWorld[slaveID].getGlobalI();
    _matA    << 0., -contact._masterR(2), contact._masterR(1),
                contact._masterR(2), 0., -contact._masterR(0),
                -contact._masterR(1), contact._masterR(0), 0.;
    _matB    << 0., -contact._slaveR(2), contact._slaveR(1),
                contact._slaveR(2), 0., -contact._slaveR(0),
                -contact._slaveR(1), contact._slaveR(0), 0.;
    _colMat  = (1./masterMass + 1./slaveMass ) * Identity
             - _matA * _masterI.inverse() * _matA - _matB * _slaveI.inverse() * _matB;
    _mass    = 1. / ( 1./masterMass + 1./slaveMass );
  }

  WrappedContact(const WrappedContact & other){
    _contact = other._contact;
    _masterI = other._masterI;
    _slaveI  = other._slaveI;
    _matA    = other._matA;
    _matB    = other._matB;
    _colMat  = other._colMat;
    _mass    = other._mass;
  }

  WrappedContact & operator=(const WrappedContact & other){
    _contact = other._contact;
    _masterI = other._masterI;
    _slaveI  = other._slaveI;
    _matA    = other._matA;
    _matB    = other._matB;
    _colMat  = other._colMat;
    _mass    = other._mass;
    return *this;
  }
};

struct PriorityInfo{
public:
  int      _contactId;
  double   _weightedV;
  double   _normalV;
  double   _energy;

  PriorityInfo(){
    _contactId = -1;
    _weightedV = 0.;
    _normalV   = 0.;
    _energy    = 0.;
  }

  PriorityInfo(int contactId, double weightedV, double normalV, double energy):
    _contactId(contactId), _weightedV(weightedV), _normalV(normalV), _energy(energy){}

  PriorityInfo(const PriorityInfo & other){
    _contactId = other._contactId;
    _weightedV = other._weightedV;
    _normalV   = other._normalV;
    _energy    = other._energy;
  }

  PriorityInfo & operator=(const PriorityInfo & other){
    _contactId = other._contactId;
    _weightedV = other._weightedV;
    _normalV   = other._normalV;
    _energy    = other._energy;
    return *this;
  }
};

struct CompareMinV{
public:
  bool operator()(const PriorityInfo & lhs, const PriorityInfo & rhs){
    return lhs._weightedV > rhs._weightedV;
    //return lhs._normalV > rhs._normalV;
  }
};

struct CompareEnergy{
public:
  bool operator()(const PriorityInfo & lhs, const PriorityInfo & rhs){
    return lhs._energy < rhs._energy;
  }
};

struct DegreeQueue{
public:
  int _id;
  int _degree;

  DegreeQueue(): _id(-1), _degree(0){}

  DegreeQueue(int id, int degree): _id(id), _degree(degree){}

  DegreeQueue(const DegreeQueue & other): _id(other._id), _degree(other._degree){}

  DegreeQueue & operator=(const DegreeQueue & other){
    _id = other._id;
    _degree = other._degree;
    return *this;
  }
};

struct CompareDegree{
public:
  bool operator()(const DegreeQueue & lhs, const DegreeQueue & rhs){
    return lhs._degree > rhs._degree;
  }
};

struct CompareSet{
public:
  bool operator()(const int & lhs, const int & rhs){
    return lhs < rhs;
  }
};

bool collisionDebug = false;
const int checkStep = 1529;
const int checkCount = 5000;
void collision_resolution(vector<ContactInfo> & contactList, int step){
  vector<WrappedContact> wrappedContactList;
  for(int i = 0; i < contactList.size(); ++i){
    if(contactList[i]._normalV < -threshold){
      int masterID = contactList[i]._master;
      int slaveID = contactList[i]._slave;
      wrappedContactList.emplace_back(contactList[i]);
    }
  }
  PriorityInfo firstItem;     //first PriorityInfo in priority_queue
  double       deltaV;        //should be positive
  double       deltaE;        //should be positive in compression->increase elastic energy
                              //but can be negative in release, because no guarantee relative
                              //velocity will be positive
  double       scalarNormalP; // impulse in normal direction, should be positive
  const int    ALPHA = 5;     // control deltaV, deltaW
  const double TOL   = 0.2;  // similarity
  const double grainMu = 0.65;  // consider friction
  const double strongeCoefficient = 0.5;
  Vector3d     impulse;       // impulse = scalarNormalP*normal + tangenP
  Vector3d     normalP;       // scalarNormalP*normal
  Vector3d     tangenP;       // impulse - scalarNormalP*normal
  Vector3d     tangDir;

  /*build priority queue*/
  std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareMinV> minVQueue;
  std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareEnergy> maxEQueue;
  /*before collision resolution, build a priority queue for min normal relative velocity*/
  for(int i = 0; i < wrappedContactList.size(); ++i){
    auto & item = wrappedContactList[i];
    if(item._contact._normalV < -threshold){ minVQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy); }
  }
  if(minVQueue.empty()){ return; }

  int count = 0;
  /*start collision phase*/
  while(!minVQueue.empty()){
    /*compression phase*/
    while(!minVQueue.empty()){
      if( collisionDebug && count > checkCount && step == checkStep){
        cout<<"############ DEBUG AT COMPRESSION PHASE: "<<count<<" #############"<<endl;
      }
      auto firstContact = wrappedContactList[minVQueue.top()._contactId];
      /*find contact with similar contact normal*/
      vector<WrappedContact> contactGroup;
      double vmin = firstContact._contact._normalV;

      int contactCheck = 1;
      while(! minVQueue.empty()){
        auto item = wrappedContactList[minVQueue.top()._contactId];
        if( item._contact._master == firstContact._contact._master ){
          contactGroup.push_back(item);
        }
        minVQueue.pop();
      }

      for(auto & contact : contactGroup){
        /*deltaV should be positive*/
        deltaV = abs( vmin ) / ALPHA / contactGroup.size();
        //if(deltaV > abs(contact._contact._normalV)){ deltaV = abs(contact._contact._normalV); }
        int slaveID = contact._contact._slave;
        int masterID = contact._contact._master;
        double slaveMass = grainsWorld[slaveID].getMass();
        double masterMass = grainsWorld[masterID].getMass();
        if(massScaling) masterMass = masterMass < slaveMass / 10. ? slaveMass / 10. : masterMass;
        /*case 1, if the friction is static*/
        impulse = deltaV * contact._colMat.inverse() * contact._contact._normal;
        scalarNormalP = impulse.dot(contact._contact._normal);
        /*normal component*/
        normalP = scalarNormalP * contact._contact._normal;
        /*tangential component*/
        tangenP = impulse - normalP;
        if(tangenP.norm() == 0){ tangDir << 0.,0.,0.; }
        else{ tangDir = tangenP/tangenP.norm(); }
        /*case 2, if the contact is dynamic*/
        if(tangenP.norm() > grainMu * normalP.norm()){
          scalarNormalP = deltaV / (contact._contact._normal.transpose() * contact._colMat * (contact._contact._normal + mu * tangDir));
          impulse = scalarNormalP * contact._contact._normal + grainMu * scalarNormalP * tangDir;
        }
        /*update velocity and omega in global frame*/
        if( collisionDebug && count > checkCount && step == checkStep){
          cout<<"contact master: "<<contact._contact._master<<" slave: "<<contact._contact._slave<<" vmin: "<<vmin<<endl;
          cout<<"master mass: "<<masterMass<<" slave mass: "<<slaveMass<<endl;
          cout<<"impulse is: "<<impulse(0)<<" "<<impulse(1)<<" "<<impulse(2)<<endl;
          Vector3d masterVel = grainsWorld[masterID].getVelocity();
          cout<<"master velocity: "<<masterVel(0)<<" "<<masterVel(1)<<" "<<masterVel(2)<<endl;
          Vector3d slaveVel = grainsWorld[slaveID].getVelocity();
          cout<<"slave velocity: "<<slaveVel(0)<<" "<<slaveVel(1)<<" "<<slaveVel(2)<<endl;
          Vector3d normaltp = contact._contact._normal;
          cout<<"weightedV: "<<contact._mass*contact._contact._normalV;
          cout<<" master: "<<contact._contact._master<<" slave: "<<contact._contact._slave<<" normal: "<<normaltp(0)<<" "<<normaltp(1)<<" "<<normaltp(2)
              <<" relative velocity: "<<contact._contact._normalV<<endl<<endl;
        }
        grainsWorld[masterID].addVelocity(impulse/masterMass);
        grainsWorld[masterID].addOmegaGlobal(contact._masterI.inverse() * contact._contact._masterR.cross(impulse));
        grainsWorld[slaveID].addVelocity(-impulse/slaveMass);
        grainsWorld[slaveID].addOmegaGlobal(contact._slaveI.inverse() * contact._contact._slaveR.cross(-impulse));

        /*update energy*/
        deltaE = 0.5 * abs(2*contact._contact._normalV + deltaV) * scalarNormalP;
        contact._contact.addEnergy(deltaE);
      }
      /*reset priority_queue*/
      std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareMinV>().swap(minVQueue);
      for(int i = 0; i < wrappedContactList.size(); ++i){
        auto & item = wrappedContactList[i];
        item._contact.computeRelVelocityThreshold();
        if(item._contact._normalV < -threshold){ minVQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy); }
      }
      count ++;
    }
    /*reset counter*/
    count = 0;
    /*restitution phase*/
    for(auto & wrappedContact : wrappedContactList){
      wrappedContact._contact._energy *= strongeCoefficient;
    }

    /*generate maxEQueue for energy release phase*/
    std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareEnergy>().swap(maxEQueue);
    for(int i = 0; i < wrappedContactList.size(); ++i){
      auto & item = wrappedContactList[i];
      item._contact.computeRelVelocityThreshold();
      if(item._contact._energy > threshold){ maxEQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy); }
    }

    /*energy release phase, much easier to converge; do not need criterion */
     while(! maxEQueue.empty() ){
       {
         auto contact = wrappedContactList[maxEQueue.top()._contactId]; maxEQueue.pop();
         double vmax = contact._contact._normalV;
         deltaV = (-vmax + sqrt(vmax*vmax +  2*contact._contact._energy * contact._contact._normal.transpose() * contact._colMat * contact._contact._normal));
         deltaV = deltaV / ALPHA;
         int slaveID = contact._contact._slave;
         int masterID = contact._contact._master;
         double slaveMass = grainsWorld[slaveID].getMass();
         double masterMass = grainsWorld[masterID].getMass();
         if(massScaling){ masterMass = masterMass < slaveMass / 10. ? slaveMass / 10. : masterMass; }
         /*consider friction*/
         /*case 1, if the friction is static*/
         impulse = deltaV * contact._colMat.inverse() * contact._contact._normal;
         scalarNormalP = impulse.dot(contact._contact._normal);
         normalP = scalarNormalP * contact._contact._normal;
         tangenP = impulse - normalP;
         if(tangenP.norm() == 0){ tangDir << 0.,0.,0.; }
         else{ tangDir = tangenP/tangenP.norm(); }
         /*case 2, if the contact is dynamic*/
         if(tangenP.norm() > grainMu * normalP.norm()){
           scalarNormalP = deltaV / (contact._contact._normal.transpose() * contact._colMat * (contact._contact._normal + mu * tangDir));
           impulse = scalarNormalP * contact._contact._normal + grainMu * scalarNormalP * tangDir;
         }
         /*update velocity and omega in global frame*/
         grainsWorld[masterID].addVelocitySeparation(impulse/masterMass);
         grainsWorld[masterID].addOmegaGlobal(contact._masterI.inverse() * contact._contact._masterR.cross(impulse));
         grainsWorld[slaveID].addVelocitySeparation( (-impulse)/slaveMass);
         grainsWorld[slaveID].addOmegaGlobal(contact._slaveI.inverse() * contact._contact._slaveR.cross(-impulse));
         /*update energy*/
         deltaE = 0.5 * (2*contact._contact._normalV + deltaV) * scalarNormalP;
         contact._contact.releaseEnergy(deltaE);
       }

      /*reset priority_queue*/
      std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareEnergy>().swap(maxEQueue);
      for(int i = 0; i < wrappedContactList.size(); ++i){
        auto & item = wrappedContactList[i];
        item._contact.computeRelVelocityThreshold();
        if(item._contact._energy > threshold){ maxEQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy); }
      }
      count ++;
    }
    count = 0;
    /*update minVQueue again, some velocities can be negative after this step*/
    std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareMinV>().swap(minVQueue);
    for(int i = 0; i < wrappedContactList.size(); ++i){
      auto & item = wrappedContactList[i];
      item._contact.computeRelVelocityThreshold();
      if(item._contact._normalV < -threshold){ minVQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy); }
    }
  }
}

#endif /* MPI_COLLISION_RESOLUTION_HPP_ */
