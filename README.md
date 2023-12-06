# 1D Finite Element Method (FEM) Solver

## Description
This project implements a one-dimensional Finite Element Method (FEM) solver in Python. It utilizes the Galerkin method with 1D Lagrange basis functions and 2nd order Gaussian Quadrature for numerical integration. The solver is capable of handling problems using Forward Euler (FE) and Backward Euler (BE) methods for time-stepping.

This project was completed for COE 352 (Advanced Scientific Computing) at the University of Texas at Austin. 

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

## Project Answers
### Question 2 - Forward Euler
Solve first by using a forward Euler time derivative discretization with a time-step of Δ𝑡 = 1/551 . Plot the results at the final time. Increase the time-step until you find the instability. What dt does this occur at? How does the solution change as N decreases?

Here are the final results, plotted at the final time, using the Forward Euler method with 551 time steps and 11 nodes.


<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/FE_11n_551t.png" alt="Your Image Description">
</p>

Here is a gif that shows many different time steps:

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/solution_evolution_0_600_free.gif" alt="Your Image Description">
</p>

Slowed down and focusing on the interval of 265 -  285 time steps, we can see that the FE method gains stability around a time step size of 1/278, or 278 time steps:

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/solution_evolution_265_285_free.gif" alt="Your Image Description">
</p>

One a fixed interval, the stability looks like this:

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/solution_evolution_265_285_fixed.gif" alt="Your Image Description">
</p>

We can also see how the number of nodes, N, affects the solution. Here is a gif that shows how the plot looks as N increases from 2 to 11.

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/solution_evolution_node_change_FE.gif" alt="Your Image Description">
</p>

As you can see, the number of nodes greatly affects the accuracy of our numerical method, with the number of nodes being directly inversly porportional to the error between the the analytical and numerical plots. 

### Question 3 -  Backward Euler
Solve the same problem with the same time-steps using an implicit backward Euler. What happens as the time-step is equal to or greater than the spatial step size? Explain why.

Here are the final results, plotted at the final time, using the Backward Euler method with 551 time steps and 11 nodes.

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/BE_11n_551t.png" alt="Your Image Description">
</p>

We can also see how the number of nodes, N, affects the solution. Here is a gif that shows how the plot looks as N increases from 2 to 11.

<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/solution_evolution_node_change_BE.gif" alt="Your Image Description">
</p>

Lastly, here is a picture of a plot when the time-step is greater than the spatial step size
<p align="center">
  <img src="https://github.com/jthet/FEM1D/blob/main/resources/BE_smallSpace_bigTime.png" alt="Your Image Description">
</p>


When the time-step Δt is large, specifically when it's comparable to or larger than the spatial step size, the solution's accuracy can degrade. This happens because the backward Euler method, while stable, is only first-order accurate in time, and large time-steps can lead to significant temporal discretization errors. As Δt increases, you we can observe that that the solution becomes less responsive to rapid changes in the solution, which is a characteristic of the method's low temporal accuracy. 
