import numpy as np
import matplotlib.pyplot as plt

'''
Jackson Thetford
Final Project for COE 352: Advanced Scientific Computing @ The University of Texas at Austin
'''

def f(x, t):
    """
    Function to represent the original problem function.

    Args:
    x (float): Spatial variable.
    t (float): Time variable.

    Returns:
    float: Value of the function at (x, t).
    """
    return ((np.pi**2 - 1)*np.exp(-t)*np.sin(np.pi*x))

def u_analyt(x, t):
    """
    Analytical solution of the differential equation.

    Args:
    x (float): Spatial variable.
    t (float): Time variable.

    Returns:
    float: Analytical solution at (x, t).
    """
    return (np.exp(-t) * np.sin(np.pi * x))

def u_init(x):
    """
    Initial condition of the differential equation.

    Args:
    x (float): Spatial variable.

    Returns:
    float: Initial condition value at x.
    """
    return np.sin(np.pi * x)

def get_user_input():
    """
    Prompt user to input the number of nodes, timesteps, and the method (FE/BE).

    Returns:
    tuple: Number of nodes (int), number of timesteps (int), method (str).
    """
    while True:
        try:
            N = int(input("How many nodes?\n"))
            nt = int(input("How many timesteps? t = [0, 1]\n"))
            method = input("Type 'FE' or 'BE' for forward/backward Euler solver methods:\n").upper()
            if N > 0 and nt > 0 and method in ['FE', 'BE']:
                break
            else:
                print("Invalid input. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a valid integer.")
    return N, nt, method

def create_mesh(N, xl=0, xr=1):
    """
    Create a computational mesh.

    Args:
    N (int): Number of nodes.
    xl (float): Left boundary. Default is 0.
    xr (float): Right boundary. Default is 1.

    Returns:
    tuple: Mesh step size (h), node locations (x), element to node map (iee).
    """
    h = (xr - xl) / (N - 1)
    x = np.linspace(xl, xr, N)
    iee = np.vstack((np.arange(N - 1), np.arange(1, N))).T
    return h, x, iee

def basis_functions(h):
    """
    Generate basis functions and their derivatives for finite element method.

    Args:
    h (float): Mesh step size.

    Returns:
    tuple: List of basis functions (phi), dx/dz, dz/dx, derivatives of phi (dphi).
    """
    phi = [0,0]
    phi[0] = lambda z: (1-z)/2
    phi[1] = lambda z: (1+z)/2
    dxdz = h/2
    dzdx = 2/h
    dphi = [-1/2, 1/2]
    return phi, dxdz, dzdx, dphi

def get_matricies(N, nnElem, nt, h, iee):
    """
    Generate global stiffness, mass, and force matrices.

    Args:
    N (int): Number of nodes.
    nnElem (int): Number of nodes per element.
    nt (int): Number of timesteps.
    h (float): Mesh step size.
    iee (array): Element to node mapping.

    Returns:
    tuple: Stiffness matrix (K), mass matrix (M), force matrix (F).
    """
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    F = np.zeros((N, nt + 1))

    klocal = np.zeros((nnElem, nnElem))
    mlocal = np.zeros((nnElem, nnElem))

    phi, dxdz, dzdx, dphi  = basis_functions(h)
    vals = [-1.0 / (3**0.5), 1.0 / (3**0.5)] # quad values
    weights = [1, 1] # quad weights

    ts = np.linspace(0, 1, nt+1)   

    for k in range(N-1):
        F[k, :] = (-1/8) * (f((vals[1]), ts) * phi[0](vals[1]) + f((vals[0]), ts) * phi[0](vals[0]))

        for l in range(nnElem): # 2 nodes per element, nnElem = number of nodes per element = 2
            global_node1 = int(iee[k][l])
            for m in range(nnElem):
                klocal[l][m] = (dphi[l]*dzdx*dphi[m]*dzdx*dxdz) * 2 # Quadrature but not dependant of time
                mlocal[l][m] = (weights[0]*phi[l](vals[0])*phi[m](vals[0]) + weights[1]*phi[l](vals[1])*phi[m](vals[1])) * dxdz * 2 # quadrature
                global_node2 = int(iee[k][m])
                K[global_node1][global_node2] += klocal[l][m]
                M[global_node1][global_node2] += mlocal[l][m]
    return K, M, F

def initialize_system( N, nnElem, nt, h, iee):
    """
    Initialize system matrices.

    Args:
    N (int): Number of nodes.
    nnElem (int): Number of nodes per element.
    nt (int): Number of timesteps.
    iee (array): Element to node mapping.

    Returns:
    tuple: Stiffness matrix (K), mass matrix (M), force matrix (F).
    """
    K, M, F = get_matricies(N, nnElem, nt, h, iee)

    # Boundary Conditions <- Reducing system
    u_db = 0 # U = 0 at x = 0, 1

    for i in range(N):
        if i == 0 or i == N-1:
            for j in range(N):
                if i != j:
                    M[j][i] = 0 
                    M[i][j] = 0
                    F[j] = F[j] - K[j][i]*u_db
            M[i][i] = 1
    return K, M, F 

def get_dirchlet(N):
    """
    Generate Dirichlet boundary condition matrix.
    u(0, t) = u(1, t) = 0

    Args:
    N (int): Number of nodes.

    Returns:
    ndarray: Dirichlet boundary condition matrix.
    """
    dbc = np.eye(N)
    dbc[0,0] = 0
    dbc[N-1,N-1] = 0
    return dbc

def solve_system(N, nt, x, M, K, F, iee, method):
    """
    Solve the system using Euler methods.

    Args:
    N (int): Number of nodes.
    nt (int): Number of timesteps.
    x (array): Node locations.
    M (ndarray): Mass matrix.
    K (ndarray): Stiffness matrix.
    F (ndarray): Force matrix.
    iee (array): Element to node mapping.
    method (str): Solver method ('FE' or 'BE').

    Returns:
    ndarray: Solution matrix.
    """
    dt = 1/nt
    invM = np.linalg.inv(M)
    M_invK = np.dot(invM, K)
    B = (1/dt) * M + K
    invB = np.linalg.inv(B)

    u = np.zeros((N, nt + 1))
    u[:, 0] = u_init(x)

    # Making dirchlet boundary conditions
    dbc = get_dirchlet(N)

    if method == 'FE':
        for t in range(nt):
            u[:, t + 1] = dbc.dot(u[:, t]-dt*M_invK.dot(u[:, t])+dt*invM.dot(F[:, t]))
    else:
        for t in range(nt):
            u[:, t + 1] = dbc.dot((1/dt)*invB.dot(M.dot(u[:, t]))+invB.dot(F[:, t]))
    return u

def line():
    """
    Print a line for separation.
    """
    print('-'*50)

def plot(u, N, timestep, nt, method):
    """
    Plot the analytical and numerical solutions.

    Args:
    u (ndarray): Solution matrix.
    N (int): Number of nodes.
    timestep (int): Timestep to plot.
    method (str): Solver method used.
    """
    x_user = np.linspace(0, 1, N)
    x_analyt = np.linspace(0, 1, 1000)
    plt.plot(x_analyt, u_analyt(x_analyt, timestep/nt), label='Analytical Solution')
    plt.plot(x_user, u[:, timestep], label=f'User Solution ({method})')
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    plt.title(f'Comparison of Analytical and Numerical Solutions using {method}')
    plt.legend()
    info_text = f'Nodes: {N}\nTime Step: {timestep}\nTime Range: [0, 1]'
    plt.text(0.05, 0.95, info_text, transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')
    plt.show()

def welcome():
    """
    Print a welcome message.
    """
    print("Welcome to the 1D FEM Solver!")
    print("\nThis code uses the Galerkin method with 1D Lagrange basis functions and 2nd order Gaussian Quadrature")

def main():
    nnElem = 2 # Nodes per element
    t = 1 # initial time
    xl = 0
    xr = 1    

    line(); welcome()
    line(); print("Set up your system . . ."); line()

    N, nt, method = get_user_input()
    h, x, iee = create_mesh(N, xl, xr)
    K, M, F = initialize_system( N, nnElem, nt, h, iee)
    u = solve_system(N, nt, x, M, K, F, iee, method)

    timestep_to_plot = nt # using the last time step for plot
    print(f'You have selected a {N} node system with {nt} time steps')
    print(f'The system is being solved with the {method} method and plotted at timestep {timestep_to_plot}')
    
    # print('K', K)
    # print('M', M)
    # print('F', F)
    # print('u', u)
    
    plot(u, N, timestep_to_plot, nt, method)

if __name__ == "__main__":
    main()
