/*
 * MPI_SOLVER_H_
 * Created on: June 21, 2020
 *      Author: Peng TAN (Berkeley)
 */
#ifndef MPI_SOLVER_H_
#define MPI_SOLVER_H_

const int T = 1;

void simulate_one_step(double dt, int step, string indir){
  stressVoigt << 0.,0.,0.,0.,0.,0.;

  for(int i = 0; i < num_grains; ++i){
    for(int j = i+1; j < num_grains; ++j){
      if(grainsWorld[i].bCircleCheck(grainsWorld[j])){
        grainsWorld[i].findInterparticleForceMoment(grainsWorld[j],dt,stressVoigt);
      }
    }
  }

  for(int i = 0; i < num_grains; ++i){

    if(step % T == 0){
      wallBalls->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
      externalWall->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
    }
    /*
    if(wallCylinder->bCircleCheck(grainsWorld[i])){
      wallCylinder->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
    }
    */
    if(wallCap->bCircleCheck(grainsWorld[i])){
      wallCap->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
    }
    if(wallPlane->bCircleCheck(grainsWorld[i])){
      wallPlane->findWallForceMoment(grainsWorld[i],dt,stressVoigt);
    }
  }

  nowCapPressure = abs( wallCap->getForce().dot(wallCap->getNormalGlobal()) )/ wallCap->getArea();
  planePressure  = abs( wallPlane->getForce().dot(wallPlane->getNormalGlobal()) )/ wallCap->getArea();
  cout<<"current cap pressure is: "<<nowCapPressure/scalingFactor<<" last cap pressure is: "<<prevCapPressure/scalingFactor<<" planePressure is: "<<planePressure/scalingFactor<<endl;
  prevCapPressure = nowCapPressure;

  for(int i = 0; i < num_grains; ++i){
    grainsWorld[i].takeTimeStep(gDamping, dt);
  }

  if(step < tConsolidate){
    wallCap->takeTimeStepConsolidate(pressure, dt);
  }else{
    wallCap->takeTimeStepStrainControl(Vector3d(0.,0.,-dHeight),dt);
  }
  wallPlane->takeTimeStepStrainControl(Vector3d(0.,0.,0.),dt);

  wallBalls->rotateTopLayer(wallCap);
  wallBalls->takeTimeStep(dt, "internalWall");
  //externalWall->rotateTopLayer(wallCap);
  //externalWall->takeTimeStep(dt, "externalWall");
  if(step >= tConsolidate) {
      wallBalls->moveHeight(-dHeight);
      //externalWall->moveHeight(-dHeight);
  }

  double kinematicEnergy = 0.0;
  double maxVelocity = 0.0;
  for(int i = 0; i < num_grains; ++i){
      Vector3d vel = grainsWorld[i].getVelocity();
      maxVelocity = maxVelocity > vel.norm() ? maxVelocity : vel.norm();
      kinematicEnergy += grainsWorld[i].computeKineticEnergy();
  }

  stressVoigt = stressVoigt/wallBalls->findVolume()/scalingFactor;
  cout<<"current stress is: "<<stressVoigt(0)<<" "<<stressVoigt(1)<<" "<<stressVoigt(2)<<" "<<stressVoigt(3)<<" "<<stressVoigt(4)<<" "<<stressVoigt(5)<<endl;
  Vector3d capNormal = wallCap->getNormalGlobal();
  cout<<"Max velocity is: "<<maxVelocity<<" Cap normal is: "<<capNormal(0)<<" "<<capNormal(1)<<" "<<capNormal(2)<<endl;

  if(step % 10 == 0){
    cout<<"THIS IS ITERATION: "<<step<<" kinematic energy is: "<<kinematicEnergy*scalingFactor*scalingFactor<<endl;
    cout<<"Volume is: "<<wallBalls->findVolume()<<" stress_33 is: "<<stressVoigt(2)<<" stress_11 is: "<<stressVoigt(0)<<endl;
    cout<<"############## END "<<step<<" ITERATIONS ################"<<endl<<endl;
  }

}
#endif /*MPI_SOLVER_H_*/
