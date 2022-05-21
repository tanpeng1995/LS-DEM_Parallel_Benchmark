#ifndef MPI_QUICK_HULL_HPP_
#define MPI_QUICK_HULL_HPP_
#include "definitions.hpp"

std::vector<std::pair<int, Vector2d>> find_hull(
  std::vector<std::pair<int, Vector2d>> & S, std::pair<int, Vector2d> P, std::pair<int, Vector2d> Q);

std::vector<std::pair<int, Vector3d>> quick_hull(std::vector<std::pair<int, Vector3d>> & S, int rank){
  /*project to a plane*/
  /*mean position of Vector2d*/
  /*
  Vector2d mean = Vector2d(0.,0.);
  for(auto item : S){
    mean += Vector2d(item.second(0), item.second(1));
  }
  mean /= S.size();
  Eigen::MatrixXd dataMat(S.size(), 2);
  for(int i = 0; i < S.size(); ++i){
    dataMat.row(i) = Vector2d( S[i].second(0), S[i].second(1) ) - mean;
  }
  Eigen::JacobiSVD<Matrix2d> svd(dataMat.transpose() * dataMat, Eigen::ComputeFullU | Eigen::ComputeFullU );
  Vector2d newX = svd.matrixU().col(svd.matrixU().cols()-1);
  Vector2d newY = Vector2d( -newX(1), newX(0) );
  newY = newY / newY.norm();
  //cout<<"newX is: "<<newX(0)<<" "<<newX(1)<<" newY is: "<<newY(0)<<" "<<newY(1)<<endl;
  std::vector<std::pair<int,Vector2d>> newS; //projected S
  for(int i = 0; i < S.size(); ++i){
    Vector2d pt = Vector2d(S[i].second(0), S[i].second(1));
    newS.push_back(std::make_pair(i, Vector2d(pt.dot(newX), pt.dot(newY))));
  }
  */
  std::vector<std::pair<int,Vector2d>> newS; //projected S
  for(int i = 0; i < S.size(); ++i){
    Vector2d pt = Vector2d(S[i].second(0), S[i].second(1));
    newS.push_back(std::make_pair(i, pt));
  }
  /*find the maximum and minimum in new projected direction*/
  std::pair<int, Vector2d> A, B;
  double leftX = std::numeric_limits<double>::max();
  double rightX = -std::numeric_limits<double>::max();
  for(auto pair : newS){
    if(pair.second(0) < leftX){
      leftX = pair.second(0);
      A = pair;
    }
    if(pair.second(0) > rightX){
      rightX = pair.second(0);
      B = pair;
    }
  }
  //cout<<"leftX is: "<<leftX<<" rightX is: "<<rightX<<endl;
  std::vector<std::pair<int, Vector2d>> convexHull;
  convexHull.push_back(A);
  convexHull.push_back(B);
  /*find normal of A and B*/
  Vector2d normal;
  normal << B.second(0) - A.second(0), A.second(1) - B.second(1);
  normal = normal / normal.norm();
  /*find S1 and S2, S1 > 0, S2 < 0*/
  std::vector<std::pair<int, Vector2d>> S1, S2;
  Vector2d diff;
  double dist;
  for(auto pair : newS){
    if(pair == A || pair == B ){ continue; }
    diff << pair.second(0)-A.second(0), pair.second(1)-A.second(1);
    dist = diff.dot(normal);
    if(dist > 0){ S1.push_back(pair); }
    else if(dist < 0){ S2.push_back(pair); }
  }
  auto subHullS1 = find_hull(S1, A, B);
  convexHull.insert(convexHull.end(), subHullS1.begin(), subHullS1.end());
  auto subHullS2 = find_hull(S2, B, A);
  convexHull.insert(convexHull.end(), subHullS2.begin(), subHullS2.end());
  std::vector<std::pair<int, Vector3d>> ret;
  for(auto item : convexHull){
    if(item.second.norm() < 125.) continue;
    ret.push_back(S[item.first]);
  }
  return ret;
}

std::vector<std::pair<int, Vector2d>> find_hull(
  std::vector<std::pair<int, Vector2d>> & S, std::pair<int, Vector2d> P, std::pair<int, Vector2d> Q){
  std::vector<std::pair<int, Vector2d>> subHull;
  if(S.size() == 0) return subHull;
  /*find the farthest point from S*/
  Vector2d normal;
  normal << Q.second(0) - P.second(0), P.second(1) - Q.second(1);
  normal = normal / normal.norm();
  Vector2d diff;
  std::pair<int, Vector2d> C;
  double maxDist = -std::numeric_limits<double>::max(), dist;
  for(auto pair : S){
    diff << pair.second(0)-P.second(0), pair.second(1)-P.second(1);
    dist = diff.dot(normal);
    if(dist < 0){ normal = -normal; dist = -dist;}
    if(dist > maxDist){
      maxDist = dist;
      C = pair;
    }
  }
  subHull.push_back(C);
  /*find S0, S1, S2, S0 inside triangle PQC, S1 on the left of PC, S2 on the right of CQ*/
  std::vector<std::pair<int, Vector2d>> S1, S2;
  Vector2d normalPC;
  normalPC << P.second(0) - C.second(0), C.second(1) - P.second(1);
  normalPC = normalPC / normalPC.norm();
  /*check with point Q, Q should be in the opposite direction of orientedd PC*/
  diff << Q.second(0) - C.second(0), Q.second(1) - C.second(0);
  if(diff.dot(normalPC) > 0){ normalPC = -normalPC; }
  Vector2d normalCQ;
  normalCQ << Q.second(0) - C.second(0), C.second(1) - Q.second(1);
  normalCQ = normalCQ / normalCQ.norm();
  /*check with point P, P should be in the opposite direction of orientedd CQ*/
  diff << P.second(0) - C.second(0), P.second(1) - C.second(0);
  if(diff.dot(normalCQ) > 0){ normalCQ = -normalCQ; }
  /*assign S1 and S2*/
  for(auto pair : S){
    if(pair == C){ continue; }
    diff << pair.second(0)-C.second(0), pair.second(1)-C.second(1);
    if(diff.dot(normalPC) > 0) { S1.push_back(pair); }
    else if(diff.dot(normalCQ) > 0) { S2.push_back(pair); }
  }
  auto ret1 = find_hull(S1, P, C);
  subHull.insert(subHull.end(), ret1.begin(), ret1.end());
  auto ret2 = find_hull(S2, C, Q);
  subHull.insert(subHull.end(), ret2.begin(), ret2.end());
  return subHull;
}

#endif /*MPI_QUICK_HULL_HPP_*/
