ID5130 Project

Mesh Generation for a domain with arbitrary local cell size requirement (2D)

The file 'Final.cpp' contains the final version of the code.
The program generates a consistent mesh of triangular elements.

Dependencies:
	The program 'gmsh' must be installed to view the generated meshes.
	It can be downloaded from: https://gmsh.info/

Usage:
	All modifiable parameters are present near the top of the file.
	The macros 'MinDepth' and 'MaxDepth' control the minimum and maximum number of divisions made to cells.
	The macro 'NumThreads' controls the number of threads used to paralellise with OpenMP.
	The function 'MeshSize' is used for specifying mesh size requirement in the domain.
		Domain is 0 <= x <= 1 and 0 <= y <= 1

	Compile and run with:
	g++ Final.cpp -fopenmp -lm
	./a.out > O.msh
	
	Open the O.msh file using gmsh
	
	The time taken by the program is written at the end of the .msh file.
		Execution time refers to total time taken, including file I/O.
		Computation time refers to time taken for computation of cells and their connectivity. Does not include file I/O.

Examples:
	The folder named 'Examples' contains the parameters, MinDepth, MaxDepth and the function MeshSize used to generate the sample meshes mentioned in the report.


