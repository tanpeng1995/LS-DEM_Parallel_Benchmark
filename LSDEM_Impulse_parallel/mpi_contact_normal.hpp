#ifndef MPI_CONTACT_NORMAL_HPP_
#define MPI_CONTACT_NORMAL_HPP_
#include <algorithm>
#include <map>
#include "Grain3d.hpp"

struct Vertex;
struct Box{
public:
  int  _bid;
  bool _inside   = false; //if has any node on the grain surface in the box
  bool _ghost    = false; //if it is the ghost box
  Vector3d _size = Vector3d(0.,0.,0.); // box size
  std::vector< std::pair<int, Vector3d> > _contactList; //contactList of grain and node id.
  Box();
};

struct Vertex{
public:
  int  _vid;
  bool _inside   = false; // if it is inside the "cylinder"
  std::vector< std::shared_ptr<Box> > _boxes; //list of boxes
  Vertex();
};

Box::Box(){
  _bid      = -1;
  _inside   = false;
  _ghost    = false;
  _size     << 0., 0., 0.;
}

Vertex::Vertex(){
  _vid      = -1;
  _inside   = false;
  _boxes.resize(8);
  std::for_each(std::begin(_boxes), std::end(_boxes),
     [](std::shared_ptr<Box> & ptr){ ptr = std::make_shared<Box>(); });
}

std::vector<std::pair<int, Vector3d>> recursive_membrane_contact_normal(
  const std::shared_ptr<Box> & domain, Vector3d bBox, int dimX, int dimY, int dimZ);

std::vector<std::pair<int, Vector3d>> recursive_membrane_contact_normal_2(
  const std::shared_ptr<Box> & domain, Vector3d bBox, int dimX, int dimY, int dimZ);

std::vector<std::pair<int, Vector3d>> membrane_contact_normal(
  Vector3d bBox, const Vector3d & domainSize, int dimX, int dimY, int dimZ){
  auto box_id = [dimX, dimY, dimZ](int i, int j, int k){
    if(i < 0 || i >= dimX+2){ i = i < 0 ? 0 : dimX+1; }
    if(j < 0 || j >= dimY+2){ j = j < 0 ? 0 : dimY+1; }
    if(k < 0 || k >= dimZ+2){ k = k < 0 ? 0 : dimZ+1; }
    return i + j*(dimX+2) + k*(dimX+2)*(dimY+2);
  };
  auto vertex_id = [dimX, dimY, dimZ](int i, int j, int k){ return i + j*(dimX+1) + k*(dimX+1)*(dimY+1); };
  double dx = domainSize(0)/dimX;
  double dy = domainSize(1)/dimY;
  double dz = domainSize(2)/dimZ;
  bBox = bBox - Vector3d(dx, dy, dz);
  /*construct boxes*/
  std::vector<std::shared_ptr<Box>> boxes((dimX+2) * (dimY+2) * (dimZ+2)); //x,y,z, include ghost boxes
  std::for_each(std::begin(boxes), std::end(boxes),
       [](std::shared_ptr<Box> & ptr){ ptr = std::make_shared<Box>(); });
  /*construct vertex*/
  std::vector<std::shared_ptr<Vertex>> vertices((dimX+1) * (dimY+1) * (dimZ+1));
  std::for_each(std::begin(vertices), std::end(vertices),
       [](std::shared_ptr<Vertex> & ptr){ ptr = std::make_shared<Vertex>(); });

  /*assign boxes to vertices*/
  for(int vid = 0; vid < vertices.size(); ++vid){
    vertices[vid]->_vid = vid;
    int k = vid/((dimX+1)*(dimY+1));                //dimension Z
    int j = vid%((dimX+1)*(dimY+1)) / (dimX+1);     //dimension Y
    int i = vid%((dimX+1)*(dimY+1)) % (dimX+1);     //dimension X

    vertices[vid]->_boxes[0] = boxes[box_id(i,j,k)];            //0
    vertices[vid]->_boxes[1] = boxes[box_id(i+1,j,k)];          //1
    vertices[vid]->_boxes[2] = boxes[box_id(i,j+1,k)];          //2
    vertices[vid]->_boxes[3] = boxes[box_id(i+1,j+1,k)];        //3

    vertices[vid]->_boxes[4] = boxes[box_id(i,j,k+1)];          //4
    vertices[vid]->_boxes[5] = boxes[box_id(i+1,j,k+1)];        //5
    vertices[vid]->_boxes[6] = boxes[box_id(i,j+1,k+1)];        //6
    vertices[vid]->_boxes[7] = boxes[box_id(i+1,j+1,k+1)];      //7
  }
  /*assign vertices to boxes*/
  for(int bid = 0; bid < boxes.size(); ++bid){
    boxes[bid]->_bid = bid;
    boxes[bid]->_size << dx, dy, dz;
    int k = bid/((dimX+2) * (dimY+2));                   //dimension Z
    int j = bid%((dimX+2) * (dimY+2)) / (dimX+2);        //dimension Y
    int i = bid%((dimX+2) * (dimY+2)) % (dimX+2);        //dimension X
    if(i < 1 || i > dimX || j < 1 || j > dimY || k < 1 || k > dimZ){
      boxes[bid]->_ghost = true;
    }
  }
  /*add node into boxes*/
  int _idx, _idy, _idz;
  int _bid, _vid;
  Vector3d pos;
  std::vector<Vector3d> pointList;
  /*normal grain*/
  for(int gid = 0; gid < num_grains; ++gid){
    if(belongToThisRank[gid] == 1){
      /*center of grain*/
      pos = grainsWorld[gid].getPosition();
      _idx = (pos(0)-bBox(0))/dx;
      _idy = (pos(1)-bBox(1))/dy;
      _idz = (pos(2)-bBox(2))/dz;
      _bid = box_id(_idx, _idy, _idz);
      if(! boxes[_bid]->_inside ) boxes[_bid]->_inside = true;
      /*surface nodes of grain*/
      pointList = grainsWorld[gid].getPointList();
      for(int i = 0; i < pointList.size(); ++i){
        pos  = pointList[i];
        _idx = (pos(0)-bBox(0))/dx;
        _idy = (pos(1)-bBox(1))/dy;
        _idz = (pos(2)-bBox(2))/dz;
        _bid = box_id(_idx, _idy, _idz);
        boxes[_bid]->_contactList.push_back(std::make_pair(gid, pos));
        if(! boxes[_bid]->_inside ) boxes[_bid]->_inside = true;
      }
    }
  }
  /*ghost grain*/
  for(int gid = 0; gid < num_grains; ++gid){
    if(belongToThisRank[gid] == 2){
      /*center of grain*/
      pos = grainsWorld[gid].getPosition();
      _idx = (pos(0)-bBox(0))/dx;
      _idy = (pos(1)-bBox(1))/dy;
      _idz = (pos(2)-bBox(2))/dz;
      _bid = box_id(_idx, _idy, _idz);
      if(! boxes[_bid]->_inside ) boxes[_bid]->_inside = true;
      /*surface nodes of grain*/
      pointList = grainsWorld[gid].getPointList();
      for(int i = 0; i < pointList.size(); ++i){
        pos  = pointList[i];
        _idx = (pos(0)-bBox(0))/dx;
        _idy = (pos(1)-bBox(1))/dy;
        _idz = (pos(2)-bBox(2))/dz;
        _bid = box_id(_idx, _idy, _idz);
        if(! boxes[_bid]->_inside) boxes[_bid]->_inside = true;
      }
    }
  }

  std::vector<std::shared_ptr<Box>> recursiveBoxes;
  std::vector<Vector3d> recursiveBoundings;

  /*traverse vertices to check if vertices are inside*/
  for(auto & vertex : vertices){
    bool inside = true;
    for(auto box : vertex->_boxes){
      if(box->_inside == false){
        inside = false;
        break;
      }
    }
    vertex->_inside = inside;
  }
  /*traverse boxes to find the boundary*/
  for(auto box : boxes){
    if(! box->_ghost){
      int bid = box->_bid;
      int k = bid/((dimX+2) * (dimY+2));                   //dimension Z
      int j = bid%((dimX+2) * (dimY+2)) / (dimX+2);        //dimension Y
      int i = bid%((dimX+2) * (dimY+2)) % (dimX+2);        //dimension X
      int outsideVertices = 0;
      if(vertices[vertex_id(i-1,j-1,k-1)]->_inside) outsideVertices++;    //0
      if(vertices[vertex_id(i,j-1,k-1)]->_inside) outsideVertices++;      //1
      if(vertices[vertex_id(i-1,j,k-1)]->_inside) outsideVertices++;      //2
      if(vertices[vertex_id(i,j,k-1)]->_inside) outsideVertices++;        //3
      if(vertices[vertex_id(i-1,j-1,k)]->_inside) outsideVertices++;      //4
      if(vertices[vertex_id(i,j-1,k)]->_inside) outsideVertices++;        //5
      if(vertices[vertex_id(i-1,j,k)]->_inside) outsideVertices++;        //6
      if(vertices[vertex_id(i,j,k)]->_inside) outsideVertices++;          //7

      if(outsideVertices > 0 && outsideVertices < 8){
        recursiveBoxes.push_back(box);
        recursiveBoundings.push_back( bBox + Vector3d(dx*i, dy*j, dz*k) );
      }
    }
  }

  std::vector<std::pair<int, Vector3d>> ret;
  /*return list of contacts*/
  for(int i = 0; i < recursiveBoxes.size(); ++i){
    auto r = recursive_membrane_contact_normal_2(recursiveBoxes[i], recursiveBoundings[i], 2, 2, 2);
    ret.insert(ret.end(), r.begin(), r.end());
  }
  return ret;
}


std::vector<std::pair<int, Vector3d>> recursive_membrane_contact_normal(
  const std::shared_ptr<Box> & domain, Vector3d bBox, int dimX, int dimY, int dimZ){
  /*stop criterion*/
  Vector3d domainSize = domain->_size;
  double dx = domainSize(0)/dimX;
  double dy = domainSize(1)/dimY;
  double dz = domainSize(2)/dimZ;
  bBox = bBox - Vector3d(dx, dy, dz);

  if(domainSize(0) <= 1. || domainSize(1) <= 1. || domainSize(2) <= 1. || domain->_contactList.size() <= 10){
    std::vector< std::pair<int, Vector3d> > returnContacts;
    std::map<int, int> groupCounts;
    std::map<int, Vector3d> groupContacts;
    for(auto contact : domain->_contactList){
      auto search = groupContacts.find(contact.first);
      if(search != groupContacts.end()){
        groupContacts[contact.first] += contact.second;
        groupCounts[contact.first]   += 1;
      }else{
        groupContacts[contact.first] = contact.second;
        groupCounts[contact.first]   = 1;
      }
    }
    for(auto item : groupContacts){
      returnContacts.emplace_back(item.first, item.second/groupCounts[item.first]);
    }
    return returnContacts;
  }
  auto box_id = [dimX, dimY, dimZ](int i, int j, int k){
    if(i < 0 || i >= dimX+2){ i = i < 0 ? 0 : dimX+1; }
    if(j < 0 || j >= dimY+2){ j = j < 0 ? 0 : dimY+1; }
    if(k < 0 || k >= dimZ+2){ k = k < 0 ? 0 : dimZ+1; }
    return i + j*(dimX+2) + k*(dimX+2)*(dimY+2);
  };
  auto vertex_id = [dimX, dimY, dimZ](int i, int j, int k){ return i + j*(dimX+1) + k*(dimX+1)*(dimY+1); };
  /*construct boxes*/
  std::vector<std::shared_ptr<Box>> boxes((dimX+2) * (dimY+2) * (dimZ+2)); //x,y,z, include ghost boxes
  std::for_each(std::begin(boxes), std::end(boxes),
       [](std::shared_ptr<Box> & ptr){ ptr = std::make_shared<Box>(); });
  /*construct vertex*/
  std::vector<std::shared_ptr<Vertex>> vertices((dimX+1) * (dimY+1) * (dimZ+1));
  std::for_each(std::begin(vertices), std::end(vertices),
       [](std::shared_ptr<Vertex> & ptr){ ptr = std::make_shared<Vertex>(); });
  /*assign boxes to vertices*/
  for(int vid = 0; vid < vertices.size(); ++vid){
    vertices[vid]->_vid = vid;
    int k = vid/((dimX+1)*(dimY+1));               //dimension Z
    int j = vid%((dimX+1)*(dimY+1)) / (dimX+1);     //dimension Y
    int i = vid%((dimX+1)*(dimY+1)) % (dimX+1);     //dimension X

    vertices[vid]->_boxes[0] = boxes[box_id(i,j,k)];            //0
    vertices[vid]->_boxes[1] = boxes[box_id(i+1,j,k)];          //1
    vertices[vid]->_boxes[2] = boxes[box_id(i,j+1,k)];          //2
    vertices[vid]->_boxes[3] = boxes[box_id(i+1,j+1,k)];        //3

    vertices[vid]->_boxes[4] = boxes[box_id(i,j,k+1)];          //4
    vertices[vid]->_boxes[5] = boxes[box_id(i+1,j,k+1)];        //5
    vertices[vid]->_boxes[6] = boxes[box_id(i,j+1,k+1)];        //6
    vertices[vid]->_boxes[7] = boxes[box_id(i+1,j+1,k+1)];      //7

  }
  /*assign vertices to boxes*/
  for(int bid = 0; bid < boxes.size(); ++bid){
    boxes[bid]->_bid = bid;
    int k = bid/((dimX+2) * (dimY+2));                   //dimension Z
    int j = bid%((dimX+2) * (dimY+2)) / (dimX+2);        //dimension Y
    int i = bid%((dimX+2) * (dimY+2)) % (dimX+2);        //dimension X
    if(i < 1 || i > dimX || j < 1 || j > dimY || k < 1 || k > dimZ){
      boxes[bid]->_ghost = true;
    }
  }

  /*add node into boxes*/
  int _idx, _idy, _idz;
  int _bid, _vid;
  int _nid, _gid;
  Vector3d pos;
  for(int i = 0; i < domain->_contactList.size(); ++i){
    _gid = domain->_contactList[i].first;
    pos  = domain->_contactList[i].second;
    _idx = (pos(0)-bBox(0))/dx;
    _idy = (pos(1)-bBox(1))/dy;
    _idz = (pos(2)-bBox(2))/dz;
    _bid = box_id(_idx, _idy, _idz);
    boxes[_bid]->_contactList.push_back(domain->_contactList[i]);
    if(! boxes[_bid]->_inside) boxes[_bid]->_inside = true;
  }

  std::vector<std::shared_ptr<Box>> recursiveBoxes;
  std::vector<Vector3d> recursiveBoundings;

  /*traverse vertices to check if vertices are inside*/
  for(auto & vertex : vertices){
    bool inside = true;
    for(auto box : vertex->_boxes){
      if(box->_inside == false){
        inside = false;
        break;
      }
    }
    vertex->_inside = inside;
  }

  /*traverse boxes to find the boundary*/
  for(auto box : boxes){
    if(! box->_ghost){
      int bid = box->_bid;
      int k = bid/((dimX+2) * (dimY+2));                   //dimension Z
      int j = bid%((dimX+2) * (dimY+2)) / (dimX+2);        //dimension Y
      int i = bid%((dimX+2) * (dimY+2)) % (dimX+2);        //dimension X
      int outsideVertices = 0;
      if(vertices[vertex_id(i-1,j-1,k-1)]->_inside) outsideVertices++;    //0
      if(vertices[vertex_id(i,j-1,k-1)]->_inside) outsideVertices++;      //1
      if(vertices[vertex_id(i-1,j,k-1)]->_inside) outsideVertices++;      //2
      if(vertices[vertex_id(i,j,k-1)]->_inside) outsideVertices++;        //3
      if(vertices[vertex_id(i-1,j-1,k)]->_inside) outsideVertices++;      //4
      if(vertices[vertex_id(i,j-1,k)]->_inside) outsideVertices++;        //5
      if(vertices[vertex_id(i-1,j,k)]->_inside) outsideVertices++;        //6
      if(vertices[vertex_id(i,j,k)]->_inside) outsideVertices++;          //7

      if(outsideVertices > 0 && outsideVertices < 8){
        recursiveBoxes.push_back(box);
        recursiveBoundings.push_back( bBox + Vector3d(dx*i, dy*j, dz*k) );
      }
    }
  }

  std::vector<std::pair<int, Vector3d>> ret;
  /*return list of contacts*/
  for(int i = 0; i < recursiveBoxes.size(); ++i){
    auto r = recursive_membrane_contact_normal(recursiveBoxes[i], recursiveBoundings[i], 2, 2, 2);
    ret.insert(ret.end(), r.begin(), r.end());
  }
  return ret;
}

std::vector<std::pair<int, Vector3d>> recursive_membrane_contact_normal_2(
  const std::shared_ptr<Box> & domain, Vector3d bBox, int dimX, int dimY, int dimZ){
  /*stop criterion*/
  Vector3d domainSize = domain->_size;
  double dx = domainSize(0)/dimX;
  double dy = domainSize(1)/dimY;
  double dz = domainSize(2)/dimZ;
  /*######################## END RECURSIVE FUNCTION ###########################*/
  if(domainSize(0) <= 1. || domainSize(1) <= 1. || domainSize(2) <= 1. || domain->_contactList.size() <= 10){
    std::vector< std::pair<int, Vector3d> > returnContacts;
    std::map<int, int> groupCounts;
    std::map<int, Vector3d> groupContacts;
    for(auto contact : domain->_contactList){
      auto search = groupContacts.find(contact.first);
      if(search != groupContacts.end()){
        groupContacts[contact.first] += contact.second;
        groupCounts[contact.first]   += 1;
      }else{
        groupContacts[contact.first] = contact.second;
        groupCounts[contact.first]   = 1;
      }
    }
    for(auto item : groupContacts){
      returnContacts.emplace_back(item.first, item.second/groupCounts[item.first]);
    }
    return returnContacts;
  }
  /*###################### CONTINUE RECURSIVE FUNCTION #########################*/
  auto box_id = [dimX, dimY, dimZ](int i, int j, int k){ return i + j*dimX + k*dimX*dimY; };
  /*construct boxes*/
  std::vector<std::shared_ptr<Box>> boxes(dimX * dimY * dimZ);
  std::for_each(std::begin(boxes), std::end(boxes),
       [](std::shared_ptr<Box> & ptr){ ptr = std::make_shared<Box>(); });
  /*add node into boxes*/
  int _idx, _idy, _idz;
  int _bid, _gid;
  Vector3d pos;
  for(int i = 0; i < domain->_contactList.size(); ++i){
    _gid = domain->_contactList[i].first;
    pos  = domain->_contactList[i].second;
    _idx = (pos(0)-bBox(0))/dx;
    _idy = (pos(1)-bBox(1))/dy;
    _idz = (pos(2)-bBox(2))/dz;
    _bid = box_id(_idx, _idy, _idz);
    boxes[_bid]->_contactList.push_back(domain->_contactList[i]);
    if(! boxes[_bid]->_inside) boxes[_bid]->_inside = true;
  }

  std::vector<std::shared_ptr<Box>> recursiveBoxes;
  std::vector<Vector3d> recursiveBoundings;

  /*traverse boxes to find the boundary*/
  for(auto box : boxes){
    if(box->_inside){
      int bid = box->_bid;
      int k = bid/(dimX * dimY);               //dimension Z
      int j = bid%(dimX * dimY) / dimX;        //dimension Y
      int i = bid%(dimX * dimY) % dimX;        //dimension X
      recursiveBoxes.push_back(box);
      recursiveBoundings.push_back( bBox + Vector3d(dx*i, dy*j, dz*k) );
    }
  }

  std::vector<std::pair<int, Vector3d>> ret;
  /*return list of contacts*/
  for(int i = 0; i < recursiveBoxes.size(); ++i){
    auto r = recursive_membrane_contact_normal_2(recursiveBoxes[i], recursiveBoundings[i], 2, 2, 2);
    ret.insert(ret.end(), r.begin(), r.end());
  }
  return ret;
}

#endif /*MPI_CONTACT_NORMAL_HPP_*/
