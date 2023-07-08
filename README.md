Numeric Solver for a Beam
=========================
The purpose of this project is to develop and implement a numeric solver for a governing equation for a beam that will represent the displacement of the beam at a certain lenght. Taking into consideration the forces and moments acting on the beam, I used 2nd order central differencing to estimate the second derivative, as well as the Gauss-Seidel method of solving a matrix to solve for displacement. Combined concepts of solids, calculus, and scientific computation to implement this solver.

Files
-----
`proj02.pdf` is the project writeup.  
`beamer.cpp` is the main method.  
`beam.hpp` is the header for the beam class.  
`beam.cpp` is the implementation for the beam class.  
`matrix.hpp` is the header for the matrix class.  
`matrix.cpp` is the implementation for the matrix class.  
`vector.hpp` is the header for the vector class.  
`vector.cpp` is teh implementation for the vector class.  
`timer.hpp` is a timer class to measure wall-clock execution time  
`init-plot.py` plots the Beam Displacement Analytic Solution  
`lot-plot.py` plots a log graph of the errors as a function of n
`runtime-plot.py` plots the Runtime Performance as a Function of n

Compiling and Running the Code
------------------------------
1. Type `g++ -o beamer beamer.cpp beam.cpp matrix.cpp vector.cpp` into the command line to compile the code.


2. Type `./beamer [n] [L] [EI] [q] <turbo>` into the command line to run the code.  

3. For example `./beamer 20 3 600 1500`

Additional Information
----------------------
For additional information about the project, see the proj02.pdf file.
