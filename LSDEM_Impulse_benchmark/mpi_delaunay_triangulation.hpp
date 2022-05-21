#ifndef MPI_DELAUNAY_TRIANGULATION_HPP_
#define MPI_DELAUNAY_TRIANGULATION_HPP_
#include "delaunator.hpp"
#include "definitions.hpp"

void applyPressure(const std::pair<int, Vector3d> & p1,
  const std::pair<int, Vector3d> & p2, const std::pair<int, Vector3d> & p3, const double & pressure, int rank);

void delaunay_triangulation(std::vector<std::pair<int, Vector3d>> & convexContacts, const double & pressure, int rank){
  /*mean position of Vector3d*/
  Vector3d mean = Vector3d(0.,0.,0.);
  for(auto contact : convexContacts){
    mean += contact.second;
  }
  mean /= convexContacts.size();
  /*2.5d triangulation, map to the best fit plane through linear square*/
  Eigen::MatrixXd dataMat(convexContacts.size(), 3);
  for(int i = 0; i < convexContacts.size(); ++i){
    dataMat.row(i) = convexContacts[i].second - mean;
  }
  /* compute X^T X, then SVD, the least eigen value correpsonding to the normal*/
  Eigen::JacobiSVD<Matrix3d> svd(dataMat.transpose() * dataMat, Eigen::ComputeFullU | Eigen::ComputeFullU );
  Vector3d planeNormal = svd.matrixU().col(svd.matrixU().cols()-1);
  /* project point onto plane, should check the order*/
  Vector3d oldX = Vector3d(1.,0.,0.);
  Vector3d newY = oldX.cross(planeNormal);
  newY = newY / newY.norm();
  Vector3d newX = planeNormal.cross(newY);
  newX = newX / newX.norm();
  /* x0, y0, x1, y1, ... */
  std::vector<double> coords;
  for(auto contact : convexContacts){
    Vector3d pos = contact.second;
    coords.push_back(pos.dot(newX));
    coords.push_back(pos.dot(newY));
  }
  /* do delaunay_triangulation */
  delaunator::Delaunator d(coords);
  // verteces are: d.triangles[i], d.triangles[i+1], d.triangles[i+2]
  for(size_t i = 0; i < d.triangles.size(); i+=3){
    applyPressure(convexContacts[d.triangles[i]], convexContacts[d.triangles[i+1]], convexContacts[d.triangles[i+2]], pressure, rank);
  }
  /* visualization for debug */
  string outdir = "/global/scratch/tanpeng/LSDEM_Modelling/Output/fabric_600_2/";
  FILE * surfaceObj  = fopen((outdir+"surface_fabric_600_"+std::to_string(rank)+".obj").c_str(), "w");
  for(auto contact : convexContacts){
    Vector3d position = contact.second;
    fprintf(surfaceObj, "v %.4f %.4f %.4f\n", position(0), position(1), position(2));
  }
  for(size_t i = 0; i < d.triangles.size(); i+=3){
    fprintf(surfaceObj, "f %d %d %d\n", d.triangles[i]+1,d.triangles[i+1]+1,d.triangles[i+2]+1);
  }
  fflush(surfaceObj);
}

void applyPressure(const std::pair<int, Vector3d> & p1,
  const std::pair<int, Vector3d> & p2, const std::pair<int, Vector3d> & p3, const double & pressure, int rank){
  Vector3d directedArea = (p2.second - p1.second).cross(p3.second - p1.second);
  Vector2d n1 = Vector2d(p1.second(0), p1.second(1));
  n1 = n1 / n1.norm();
  Vector2d n2 = Vector2d(directedArea(0), directedArea(1));
  n2 = n2 / n2.norm();
  if(n1.dot(n2) > 0) directedArea = -directedArea;
  Vector3d normal = directedArea / directedArea.norm();
  cout<<"rank: "<<rank<<" normal :"<<normal(0)<<" "<<normal(1)<<" "<<normal(2)<<endl;
  Vector3d pointForce   = directedArea * pressure / 6.;
  /*apply force on each grain, not implement here*/
}

#endif /*MPI_DELAUNAY_TRIANGULATION_HPP_*/
