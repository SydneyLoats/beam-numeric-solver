Numeric Solver for a Beam
=========================
The purpose of this project is to develop and implement a numeric solver for a governing equation for a beam that will represent the displacement of the beam at a certain lenght. Taking into consideration the forces and moments acting on the beam, I used 2nd order central differencing to estimate the second derivative, as well as the Gauss-Seidel method of solving a matrix to solve for displacement. Combined concepts of solids, calculus, and scientific computation to implement this solver.

Compiling and Running the Code
------------------------------
1. Type `g++ -o beamer beamer.cpp beam.cpp matrix.cpp vector.cpp` into the command line to compile the code.


2. Type `./beamer [n] [L] [EI] [q] <turbo>` into the command line to run the code.  

Additional Information
----------------------
For additional information about the project, see the proj02.pdf file.
