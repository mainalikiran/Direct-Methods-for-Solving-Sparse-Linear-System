# Research Project: Direct Methods for Solving Sparse Linear System (Lawrence Berkeley National Laboratory, CA)

## Abstract:
 Multiple methods are available in the literature to solve a linear system. They are generally classified in two categories, direct methods and iterative methods. In this work, we focus on solving a sparse linear system  _Ax = b_, where *A* is symmetric positive definite. In this case, forward and backward triangular solve can be used after the Cholesky decomposition of matrix _(A = LL^T)_
 has been computed. In this work, we explore the use of Block Low- Rank Compression(BLR) in the left-looking Cholesky algorithm. To this end, we investigate sorting the updaters to a certain target block during the factorization process and it's impact on rank growth.we compare the results with a regular factorization process without sorting the updates to see the growth of the rank of the target block.

## Introduction:
Large sparse matrices  appear in scientific applications in numerous fields including discretized PDEs, optimization problems, circuit design, structural dynamics and many more.  Solving  large sparse linear systems is an important issue in academic research and in industrial applications. The development of supercomputers has provided us with the opportunity to solve these critical problems . These challanges brought some new issues such as memory consumption, accuracy, and speed, which are encountered daily by numerical analyst and mathematicians. 
