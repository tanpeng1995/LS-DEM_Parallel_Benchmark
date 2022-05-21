# LS-DEM
This is Peng Tan's code for LS-DEM
For more details, go into each folder and there will be more explanations

## LSDEM_Parallel
- this is the major code used for simulate the triaxial compression test, this is the most reliable one
- it is parallelized, with an exclusive processor to handle membrane-grain interaction

## LSDEM_Membrane
- again, for the triaxial compression test
- it is parallelized, but model membrane in every processor

## LSDEM_Impulse_serial
- this is the serial code for impulse-based LS-DEM
- it is most feasible as it is serial, the user can change the environment setting, e.g., to handle arbitrarily shaped topography

## LSDEM_Impulse_parallel
- the parallel version of impulse-based LS-DEM with domain decomposition
- it is not effecient because domain decomposition cannot properly deal with the contact island
- our benchmarked result shows that it can at most use 8 processors, more than that will induce heavy communication overhead

## LSDEM_Impulse_benchmark
- almost the same code as LSDEM_Impulse_parallel, but used for a particular benchmarking case

## LSDEM_ball
- to study the effect of particle shape, this code simulate spherical particle with equivalent volume as LS avatars
- the same contact model and parameters are used.

## LSDEM_Settle
- for low-resolution reconstruction, need to make the particle to settle down

## LSDEM_Stablize
- for low-resolution reconstruction, the next step after settlement is to reduce internal stress generated inside the sample

## LSDEM_Test
- this is used to benchmark the performance of parallel code

## LSDEM_Test_ducplicate
- for very large parallel benchmark problem, we duplicated some avatars from existing ones, e.g, to expand the domain of size 200x200x200 to 400x400x400 by duplicating 8 copies of existing avatars


