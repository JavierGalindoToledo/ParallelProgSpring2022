# Parallel Programming Course

The following parallel programming technologies are considered in practice:
  * `MPI`
  * `OpenMP`
  * `TBB`
  * `std::thread`
## Task Topic
Fox algorithm: This algorithm breaks the input matrices down into blocks and completes the multiplication by passing these blocks around between the processes. The algorithm is designed to solve truely large multiplication problems where the input matrices are too large to be sent to each process.
## General project structure
  - `modules` Contains folders of each variant.
     -  `task_1\galindo_fox_algorithm` MPI version.
     -  `task_2\galindo_fox_algorithm_omp` OpenMP version.
     -  `task_1\galindo_fox_algorithm_tbb` Threading building blocks (TBB) version.
     -  `task_1\galindo_fox_algorithm_thread` Thread class version.
