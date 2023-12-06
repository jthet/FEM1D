import fem1d as myfem
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os


def plot(u, N, timestep, method, filename=None):
    x_user = np.linspace(0, 1, N)
    x_analyt = np.linspace(0, 1, 1000)
    # plt.plot(x_analyt, myfem.u_analyt(x_analyt, 1), label='Analytical Solution') ### use for dif step sizes
    plt.plot(x_analyt, myfem.u_analyt(x_analyt, timestep/nt), label='Analytical Solution') # use for through time
    plt.plot(x_user, u[:, timestep], label=f'User Solution ({method})')
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    plt.title(f'Comparison of Analytical and Numerical Solutions using {method}')
    plt.legend()
    info_text = f'Nodes: {N}\nTime Step: {timestep}\nTime Range: [0, 1]'
    # info_text = f'Nodes: {N}\nTime Step Size: {timestep}\nTime Range: [0, 1]'
    plt.text(0.05, 0.95, info_text, transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')
    
    # Perserving axis
    plt.xlim(0, 1)
    #plt.ylim(0, 1) # Use to show through time

    #plt.ylim(0, 0.5) # Use to show dif time step sizes at last time step
    
    if filename:
        plt.savefig(filename)
    plt.close()

def create_gif_dif_numNodes():
    global nt
    nt = 551
    num_of_steps = nt
    nnElem = 2
    xl = 0
    xr = 1    
    N = [i for i in range(2, 20, 1)]
    method = 'BE'
    # method = 'FE'

    # Create 'resources' directory if it doesn't exist
    if not os.path.exists('resources'):
        os.makedirs('resources')

    filenames = []
    for numNodes in N:
        nt = num_of_steps  # Assuming nt should vary with step_size
        h, x, iee = myfem.create_mesh(numNodes, xl, xr)
        K, M, F = myfem.initialize_system(numNodes, nnElem, nt, h, iee)
        u = myfem.solve_system(numNodes, nt, x, M, K, F, iee, method)

        filename = os.path.join('resources', f'plot_nodes_{numNodes}.png')
        plot(u, numNodes, nt, method, filename)
        filenames.append(filename)

    # Create a GIF
    gif_filename = os.path.join('resources', 'solution_evolution.gif')
    with imageio.get_writer(gif_filename, mode='I', fps = 2) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            os.remove(filename)  # Remove file after adding it to the GIF

    print(f"GIF created: {gif_filename}")
    return


def create_gif_dif_stepSizes():
    global nt
    num_of_steps = [i for i in range(1, 600, 10)]
    nnElem = 2
    xl = 0
    xr = 1    
    N = 11 # number of nodes
    # method = 'FE'
    method = 'BE'

    # Create 'resources' directory if it doesn't exist
    if not os.path.exists('resources'):
        os.makedirs('resources')

    filenames = []
    for steps in num_of_steps:
        nt = steps  # Assuming nt should vary with step_size
        h, x, iee = myfem.create_mesh(N, xl, xr)
        K, M, F = myfem.initialize_system(N, nnElem, nt, h, iee)
        u = myfem.solve_system(N, nt, x, M, K, F, iee, method)

        filename = os.path.join('resources', f'plot_stepsize_{steps}.png')
        plot(u, N, nt, method, filename)
        filenames.append(filename)

    # Create a GIF
    gif_filename = os.path.join('resources', 'solution_evolution.gif')
    with imageio.get_writer(gif_filename, mode='I', fps = 10) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            os.remove(filename)  # Remove file after adding it to the GIF

    print(f"GIF created: {gif_filename}")

def create_gif_through_time():
    global h, nnElem, nt
    nnElem = 2 # Nodes per element
    t = 1 # initial time
    xl = 0
    xr = 1    

    N = 11
    nt = 551
    # method = 'FE'
    method = 'BE'

    h, x, iee = myfem.create_mesh(N, xl, xr)
    K, M, F = myfem.initialize_system( N, nnElem, nt, h, iee)
    u = myfem.solve_system(N, nt, x, M, K, F, iee, method)

    timestep_to_plot = nt # using the last time step for plot

    if not os.path.exists('resources'):
        os.makedirs('resources')

    filenames = []
    for timestep in range(int((nt+1)/10)):
        filename = os.path.join('resources', f'frame_{timestep*10}.png')
        plot(u, N, timestep*10, method, filename)
        filenames.append(filename)

    # Create a GIF
    gif_filename = os.path.join('resources', 'solution_evolution.gif')
    with imageio.get_writer(gif_filename, mode='I', fps = 60) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            os.remove(filename)

    print(f"GIF created: {gif_filename}")

if __name__ == "__main__":
    create_gif_through_time()
    # create_gif_dif_stepSizes()
    # create_gif_dif_numNodes()

