#ifndef MPI_TRIANGULATION_PRESSURE_HPP_
#define MPI_TRIANGULATION_PRESSURE_HPP_
#include "Grain3d.hpp"
#include "definitions.hpp"
#include <algorithm>

void makeTriangle(const std::pair<int, Vector3d> & p1,
  const std::pair<int, Vector3d> & p2, const std::pair<int, Vector3d> & p3, const double & pressure);

void triangulation_apply_pressure(
  std::vector<std::vector<std::pair<int, Vector3d>>> & convexHullRings, const double & pressure){
  /*sort convexHullRings*/
  for(auto convexHullRing : convexHullRings){
    std::sort(std::begin(convexHullRing), std::end(convexHullRing),
      [](const std::pair<int, Vector3d> & lhs, const std::pair<int, Vector3d> & rhs){
        double temp1 = lhs.second(1)/(pow(lhs.second(0),2) + pow(lhs.second(1),2));
        double temp2 = rhs.second(1)/(pow(rhs.second(0),2) + pow(rhs.second(1),2));
        return temp1 < temp2;
      });
  }
  /*build triangle*/
  for(int k = 1; k < convexHullRings.size(); ++k){
    int sizeOfUpRing  = convexHullRings[k].size();
    int sizeOfLowRing = convexHullRings[k-1].size();
    double densityRatio = static_cast<double>(sizeOfLowRing-1) / static_cast<double>(sizeOfUpRing-1);
    int p = 0; //prev index of down ring
    int q = 0; //current index of down ring
    for(int i = 1; i < sizeOfUpRing; ++i){
      q = static_cast<int>(densityRatio*i);
      for(int temp = p; temp < q; ++temp){
        //make sure clockwise
        makeTriangle(convexHullRings[k][i-1], convexHullRings[k-1][temp+1], convexHullRings[k-1][temp], pressure);
      }
      //make sure clockwise
      makeTriangle(convexHullRings[k][i-1], convexHullRings[k][i], convexHullRings[k-1][q], pressure);
      p = q;
    }
  }
}

void makeTriangle(const std::pair<int, Vector3d> & p1,
  const std::pair<int, Vector3d> & p2, const std::pair<int, Vector3d> & p3, const double & pressure){
  Vector3d directedArea = (p2.second - p1.second).cross(p3.second - p1.second);
  Vector3d pointForce   = directedArea * pressure / 6.;
  /*apply force on each grain, not implement here*/
}

#endif /*MPI_TRIANGULATION_PRESSURE_HPP_*/
