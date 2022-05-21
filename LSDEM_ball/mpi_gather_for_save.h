/*
* MPI_GATHER_FOR_SAVE_H_
* Created on: June 13, 2020
*      Author: Peng TAN (Berkeley)
*/
#ifndef MPI_GATHER_FOR_SAVE_H_
#define MPI_GATHER_FOR_SAVE_H_

void gather_for_save(int step, string outdir, string indir, FILE* posfile, FILE* rotfile, FILE* wallfile, FILE* wallfile2, FILE * stressFile, FILE * volumeCNFile){
   //output
   for (int i = 0; i < num_grains; ++i) {
     const Vector3d position = grainsWorld[i].getPosition();
     const double radius = grainsWorld[i].getRadius();
     fprintf(posfile, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
     fprintf(rotfile, "%.4f\n", radius);
   }
   fflush(posfile);
   fflush(rotfile);

   for(int i = 0; i < wallBalls->getNumWalls(); ++i){
     const Vector3d position = wallBalls->getPosition(i);
     fprintf(wallfile, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
   }
   for(int i = 0; i < externalWall->getNumWalls(); ++i){
     const Vector3d position = externalWall->getPosition(i);
     fprintf(wallfile2, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
   }
   fprintf(stressFile, "%.4f %.4f %.4f %.4f %.4f %.4f\n", stressVoigt(0), stressVoigt(1), stressVoigt(2),stressVoigt(3),stressVoigt(4),stressVoigt(5));
   fflush(stressFile);
   fflush(wallfile);
   fflush(wallfile2);
   fprintf(volumeCNFile, "%.4f \n", wallBalls->findVolume()*scalingFactor*scalingFactor*scalingFactor);
   fflush(volumeCNFile);
}
#endif /*MPI_GATHER_FOR_SAVE_H_*/
