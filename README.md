Age-structured local diffusion models have important applications in population 
ecology and infectious disease dynamics, and their stability analysis typically reduces to eigenvalue 
problems of the associated infinite-dimensional linear operators. In this repository, pseudospectral 
methods are employed to study the computation of eigenvalue problems for such models. For 
various types of spatial boundary conditions (Dirichlet, Neumann, Robin, and periodic) and the 
integral boundary condition in age, we enforce the discretized boundary conditions as constraints 
appended to the differentiation matrix obtained from the discretization of the equations, ultimately 
yielding a matrix generalized eigenvalue problem. This approach provides a concise and unified
way to handle different types of boundary conditions.
