import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy.interpolate import griddata
import os

def main():
    steps = np.array([1/8, 1/16, 1/32, 1/64, 1/128])
    times = np.array(range(3, 9))

    fault = (np.loadtxt("faults_runge_naive.txt"))
    # fault = (np.loadtxt("faults_implicit.txt"))
    
    fault.shape = (len(steps), len(times))

    fig = plt.figure("FAULT")

    ax = fig.add_subplot(111, projection = '3d')

    xi = np.linspace(steps.min(),steps.max(),100)
    yi = np.linspace(times.min(),times.max(),100)
    # zi = griddata((steps, times), Z, (xi[None,:], yi[:,None]), method='cubic')
    # zi = griddata((steps, times), Z, (xi[None,:], yi[:,None]), method='cubic')

    x, y = np.meshgrid(times, steps)
    print (x)
    print (y)
    print(fault)


    ax.plot_surface(x, y, fault, alpha=0.9, rstride=1, cstride=1, linewidth=0.5, cmap='jet')
    # ax.scatter(x, y, fault)

    ax.set_xlabel('dt')
    ax.set_ylabel('h')
    ax.set_zlabel('fault')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: '1e-{:d}'.format(int(x))))

    ax.view_init(azim=-30)
    plt.savefig('plot.png', transparent=True)
    plt.show()

if __name__ == '__main__':
    main()