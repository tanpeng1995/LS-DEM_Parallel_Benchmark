#ifndef GATHER_FOR_SAVE_HPP_
#define GATHER_FOR_SAVE_HPP_

void gather_for_save(int step, string outdir, FILE* posfile, FILE* rotfile, FILE* wallfile1, FILE * wallfile2){
  //output
  for (int i = 0; i < num_grains; ++i) {
    const Vector3d position = grainsWorld[i]->getPosition();
    const Vector4d rotation = grainsWorld[i]->getQuat();
    fprintf(posfile, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
    fprintf(rotfile, "%.4f %.4f %.4f %.4f\n", rotation(0), rotation(1), rotation(2), rotation(3));
  }
  fflush(posfile);
  fflush(rotfile);
  /*wallBalls*/
  for(int i = 0; i < wallBalls[0]->getNumWalls(); ++i){
    const Vector3d position = wallBalls[0]->getPosition(i);
    fprintf(wallfile1, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
  }
  fflush(wallfile1);
  for(int i = 0; i < wallBalls[1]->getNumWalls(); ++i){
    const Vector3d position = wallBalls[1]->getPosition(i);
    fprintf(wallfile2, "%.4f %.4f %.4f\n", position(0), position(1), position(2));
  }
  fflush(wallfile2);
}
#endif /*GATHER_FOR_SAVE_HPP_*/
