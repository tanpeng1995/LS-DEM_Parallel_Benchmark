#ifndef COLLISION_RESOLUTION_HPP_
#define COLLISION_RESOLUTION_HPP_
#include "Grain3d.hpp"
#include "definitions.hpp"
#include "utilities.hpp"
#include "init_simulation.hpp"
#include "WallBalls.hpp"
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
    Identity << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    int masterId = contact._master;
    int slaveId = contact._slave;
    if(slaveId < 0){
      _masterI = grainsWorld[masterId]->getGlobalI();
      _matA    << 0., -contact._masterR(2), contact._masterR(1),
                  contact._masterR(2), 0., -contact._masterR(0),
                  -contact._masterR(1), contact._masterR(0), 0.;
      _colMat  = 1./grainsWorld[masterId]->getMass() * Identity - _matA * _masterI.inverse() * _matA;
      _mass    = grainsWorld[masterId]->getMass();
    }else if(slaveId >= num_grains){
      const int wallId = whichWallBall(slaveId);
      _masterI = grainsWorld[masterId]->getGlobalI();
      double ballMoin = wallBalls[wallId]->getMoin();
      double ballMass = wallBalls[wallId]->getMass();
      _slaveI  << ballMoin, 0., 0., 0., ballMoin, 0., 0., 0., ballMoin;
      _matA    << 0., -contact._masterR(2), contact._masterR(1),
                  contact._masterR(2), 0., -contact._masterR(0),
                  -contact._masterR(1), contact._masterR(0), 0.;
      _matB    << 0., -contact._slaveR(2), contact._slaveR(1),
                  contact._slaveR(2), 0., -contact._slaveR(0),
                  -contact._slaveR(1), contact._slaveR(0), 0.;
      _colMat  = (1./grainsWorld[masterId]->getMass() + 1./ballMass) * Identity
               - _matA * _masterI.inverse() * _matA - _matB * _slaveI.inverse() * _matB;
      _mass    = 1. / ( 1./grainsWorld[masterId]->getMass() + 1./ballMass );
    }else{
      _masterI = grainsWorld[masterId]->getGlobalI();
      _slaveI  = grainsWorld[slaveId]->getGlobalI();
      _matA    << 0., -contact._masterR(2), contact._masterR(1),
                  contact._masterR(2), 0., -contact._masterR(0),
                  -contact._masterR(1), contact._masterR(0), 0.;
      _matB    << 0., -contact._slaveR(2), contact._slaveR(1),
                  contact._slaveR(2), 0., -contact._slaveR(0),
                  -contact._slaveR(1), contact._slaveR(0), 0.;
      _colMat  = (1./grainsWorld[masterId]->getMass() + 1./grainsWorld[slaveId]->getMass()) * Identity
               - _matA * _masterI.inverse() * _matA - _matB * _slaveI.inverse() * _matB;
      _mass    = 1. / ( 1./grainsWorld[masterId]->getMass() + 1./grainsWorld[slaveId]->getMass() );
    }
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
const int checkStep = 519;
const int checkCount = 5000;
void collision_resolution(vector<ContactInfo> & contactList, int step){
  vector<WrappedContact> wrappedContactList;
  for(int i = 0; i < contactList.size(); ++i){
    if(contactList[i]._normalV < -threshold){
      wrappedContactList.emplace_back(contactList[i]);
    }
  }
  PriorityInfo firstItem;     //first PriorityInfo in priority_queue
  double       deltaV;        //should be positive
  double       deltaE;        //should be positive in compression->increase elastic energy
                              //but can be negative in release, because no guarantee relative
                              //velocity will be positive
  double       scalarNormalP; // impulse in normal direction, should be positive
  const int    ALPHA = 10;     // control deltaV, deltaW
  const double TOL   = 0.2;  // similarity
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
      {
        auto contact = wrappedContactList[minVQueue.top()._contactId];
        double vmin = contact._contact._normalV;
        /*deltaV should be positive*/
        deltaV = abs( vmin ) / ALPHA;
        int slaveId = contact._contact._slave;
        int masterId = contact._contact._master;
        double slaveMass, masterMass;
        if(slaveId < 0){
          slaveMass = std::numeric_limits<double>::max();
        }else if(slaveId >= num_grains){
          const int wallId = whichWallBall(slaveId);
          slaveMass = wallBalls[wallId]->getMass();
        }else{
          slaveMass = grainsWorld[slaveId]->getMass();
        }
        masterMass = grainsWorld[masterId]->getMass();
        /*case 1, if the friction is static*/
        impulse = deltaV * contact._colMat.inverse() * contact._contact._normal;
        scalarNormalP = impulse.dot(contact._contact._normal);
        /*normal component*/
        normalP = scalarNormalP * contact._contact._normal;
        /*tangential component*/
        tangenP = impulse - normalP;
        if(tangenP.norm() == 0){ tangDir << 0.,0.,0.; }
        else{ tangDir = tangenP/tangenP.norm(); }
        double grainMu = slaveId < 0 ? wallFriction : mu;
        /*case 2, if the contact is dynamic*/
        if(tangenP.norm() > grainMu * normalP.norm()){
          scalarNormalP = deltaV / (contact._contact._normal.transpose() * contact._colMat * (contact._contact._normal + mu * tangDir));
          impulse = scalarNormalP * contact._contact._normal + grainMu * scalarNormalP * tangDir;
        }
        /* debug */
        if( collisionDebug && count > checkCount && step == checkStep ){
          cout<<"############ DEBUG AT COMPRESSION PHASE: "<<count<<" #############"<<endl;
          cout<<"contact master: "<<contact._contact._master<<" slave: "<<contact._contact._slave<<" vmin: "<<vmin<<endl;
          cout<<"master mass: "<<masterMass<<" slave mass: "<<slaveMass<<endl;
          cout<<"impulse is: "<<impulse(0)<<" "<<impulse(1)<<" "<<impulse(2)<<endl;
          Vector3d masterVel = grainsWorld[masterId]->getVelocity();
          cout<<"master velocity: "<<masterVel(0)<<" "<<masterVel(1)<<" "<<masterVel(2)<<endl;
          if(slaveId >= 0 && slaveId < num_grains){
            Vector3d slaveVel = grainsWorld[slaveId]->getVelocity();
            cout<<"slave velocity: "<<slaveVel(0)<<" "<<slaveVel(1)<<" "<<slaveVel(2)<<endl;
          }
          Vector3d normaltp = contact._contact._normal;
          cout<<"weightedV: "<<contact._mass*contact._contact._normalV;
          cout<<" master: "<<contact._contact._master<<" slave: "<<contact._contact._slave<<" normal: "<<normaltp(0)<<" "<<normaltp(1)<<" "<<normaltp(2)
              <<" relative velocity: "<<contact._contact._normalV<<endl<<endl;
        }
        /*update velocity and omega in global frame*/
        grainsWorld[masterId]->addVelocity(impulse/masterMass);
        grainsWorld[masterId]->addOmegaGlobal(contact._masterI.inverse() * contact._contact._masterR.cross(impulse));
        if(slaveId >= 0){
          if(slaveId >= num_grains){
            const int wallId = whichWallBall(slaveId);
            wallBalls[wallId]->addVelocity(slaveId-wallBalls[wallId]->getStartIdx(), (-impulse)/slaveMass);
          }else{
            grainsWorld[slaveId]->addVelocity(-impulse/slaveMass);
            grainsWorld[slaveId]->addOmegaGlobal(contact._slaveI.inverse() * contact._contact._slaveR.cross(-impulse));
          }
        }
        /*update energy*/
        deltaE = 0.5 * abs(2*contact._contact._normalV + deltaV) * scalarNormalP;
        contact._contact.addEnergy(deltaE);
      }
      /*reset priority_queue*/
      std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareMinV>().swap(minVQueue);
      for(int i = 0; i < wrappedContactList.size(); ++i){
        auto & item = wrappedContactList[i];
        if(item._contact._slave < 0){
          item._contact.computeRelVelocityWithTopographyThreshold();
        }else if(item._contact._slave >= num_grains){
          const int wallId = whichWallBall(item._contact._slave);
          item._contact.computeRelVelocityWithBallsThreshold(wallBalls[wallId]->getVelocity(item._contact._slave - wallBalls[wallId]->getStartIdx()));
        }else{
          item._contact.computeRelVelocityThreshold();
        }
        if(item._contact._normalV < -threshold){
          minVQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy);
        }
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
      if(item._contact._slave < 0){
        item._contact.computeRelVelocityWithTopographyThreshold();
      }else if(item._contact._slave >= num_grains){
        const int wallId = whichWallBall(item._contact._slave);
        item._contact.computeRelVelocityWithBallsThreshold(wallBalls[wallId]->getVelocity(item._contact._slave - wallBalls[wallId]->getStartIdx()));
      }else{
        item._contact.computeRelVelocityThreshold();
      }
      if(item._contact._energy > threshold){
        maxEQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy);
      }
    }

    /*energy release phase*/
    while(! maxEQueue.empty() ){
       {
         auto contact = wrappedContactList[maxEQueue.top()._contactId]; maxEQueue.pop();
         double vmax = contact._contact._normalV;
         deltaV = (-vmax + sqrt(vmax*vmax +  2*contact._contact._energy * contact._contact._normal.transpose() * contact._colMat * contact._contact._normal));
         deltaV = deltaV / ALPHA;
         int slaveId = contact._contact._slave;
         int masterId = contact._contact._master;
         double slaveMass, masterMass;
         if(slaveId < 0){
           slaveMass = std::numeric_limits<double>::max();
         }else if(slaveId >= num_grains){
           const int wallId = whichWallBall(slaveId);
           slaveMass = wallBalls[wallId]->getMass();
         }else{
           slaveMass = grainsWorld[slaveId]->getMass();
         }
         masterMass = grainsWorld[masterId]->getMass();
         /*consider friction*/
         /*case 1, if the friction is static*/
         impulse = deltaV * contact._colMat.inverse() * contact._contact._normal;
         scalarNormalP = impulse.dot(contact._contact._normal);
         normalP = scalarNormalP * contact._contact._normal;
         tangenP = impulse - normalP;
         if(tangenP.norm() == 0){ tangDir << 0.,0.,0.; }
         else{ tangDir = tangenP/tangenP.norm(); }
         /*case 2, if the contact is dynamic*/
         double grainMu = slaveId < 0 ? wallFriction : mu;
         if(tangenP.norm() > grainMu * normalP.norm()){
           scalarNormalP = deltaV / (contact._contact._normal.transpose() * contact._colMat * (contact._contact._normal + mu * tangDir));
           impulse = scalarNormalP * contact._contact._normal + grainMu * scalarNormalP * tangDir;
         }
         /*update velocity and omega in global frame*/
         grainsWorld[masterId]->addVelocitySeparation(impulse/masterMass);
         grainsWorld[masterId]->addOmegaGlobal(contact._masterI.inverse() * contact._contact._masterR.cross(impulse));
         if(slaveId >= 0){
           if(slaveId >= num_grains){
             const int wallId = whichWallBall(slaveId);
             wallBalls[wallId]->addVelocity(slaveId-wallBalls[wallId]->getStartIdx(), (-impulse)/slaveMass);
           }else{
             grainsWorld[slaveId]->addVelocitySeparation( (-impulse)/slaveMass);
             grainsWorld[slaveId]->addOmegaGlobal(contact._slaveI.inverse() * contact._contact._slaveR.cross(-impulse));
           }
         }
         /*update energy*/
         deltaE = 0.5 * (2*contact._contact._normalV + deltaV) * scalarNormalP;
         contact._contact.releaseEnergy(deltaE);
       }

      /*reset priority_queue*/
      std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareEnergy>().swap(maxEQueue);
      for(int i = 0; i < wrappedContactList.size(); ++i){
        auto & item = wrappedContactList[i];
        if(item._contact._slave < 0){
          item._contact.computeRelVelocityWithTopographyThreshold();
        }else if(item._contact._slave >= num_grains){
          const int wallId = whichWallBall(item._contact._slave);
          item._contact.computeRelVelocityWithBallsThreshold(wallBalls[wallId]->getVelocity(item._contact._slave - wallBalls[wallId]->getStartIdx()));
        }else{
          item._contact.computeRelVelocityThreshold();
        }
        if(item._contact._energy > threshold){
          maxEQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy);
        }
      }
      count ++;
    }
    count = 0;
    /*update minVQueue again, some velocities can be negative after this step*/
    std::priority_queue<PriorityInfo, std::vector<PriorityInfo>, CompareMinV>().swap(minVQueue);
    for(int i = 0; i < wrappedContactList.size(); ++i){
      auto & item = wrappedContactList[i];
      if(item._contact._slave < 0){
        item._contact.computeRelVelocityWithTopographyThreshold();
      }else if(item._contact._slave >= num_grains){
        const int wallId = whichWallBall(item._contact._slave);
        item._contact.computeRelVelocityWithBallsThreshold(wallBalls[wallId]->getVelocity(item._contact._slave - wallBalls[wallId]->getStartIdx()));
      }else{
        item._contact.computeRelVelocityThreshold();
      }
      if(item._contact._normalV < -threshold){
        minVQueue.emplace(i, item._mass*item._contact._normalV, item._contact._normalV, item._contact._energy);
      }
    }
  }
}

#endif /* COLLISION_RESOLUTION_HPP_ */
