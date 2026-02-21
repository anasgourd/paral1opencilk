# Sparse Matrix Multiplication with OpenCilk

This C project implements parallel sparse matrix multiplication using CSR and CSC formats with OpenCilk. 

It is the parallel extension of a sequential implementation and was developed for the Parallel and Distributed Systems course at ECE AUTH (2023–2024).

An example run is included demonstrating correct execution and race-free behavior using the Cilksan race detector. The program can be run with any suitable matrix by providing the corresponding Matrix Market file
