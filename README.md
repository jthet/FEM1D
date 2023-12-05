# 1D Finite Element Method (FEM) Solver

## Description
This project implements a one-dimensional Finite Element Method (FEM) solver in Python. It utilizes the Galerkin method with 1D Lagrange basis functions and 2nd order Gaussian Quadrature for numerical integration. The solver is capable of handling problems using Forward Euler (FE) and Backward Euler (BE) methods for time-stepping.

## Features
- **Solver Methods**: Implements both Forward Euler and Backward Euler methods.
- **Basis Functions**: Uses 1D Lagrange basis functions for spatial discretization.
- **Gaussian Quadrature**: Employs 2nd order Gaussian Quadrature for numerical integration.
- **User Interaction**: Allows user input for nodes, timesteps, and solver methods.

## Requirements
- Python 3.x
- NumPy
- Matplotlib

## Usage
To run the solver, execute the `fem1d.py` script in a Python environment:
```
python3 fem1d.py
```


Follow the on-screen prompts to enter the number of nodes, the number of timesteps, and the preferred solver method (FE/BE).

## Input Parameters
- `N`: Number of spatial nodes in the mesh.
- `nt`: Number of timesteps for the temporal discretization.
- `method`: Choice of solver method - 'FE' for Forward Euler or 'BE' for Backward Euler.

## Output
The program computes and plots the numerical solution against the analytical solution at the specified timestep, providing a visual comparison.



