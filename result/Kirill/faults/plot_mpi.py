import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy.interpolate import griddata
import os

def main():
    nodes = np.array([1, 2, 4, 8, 16])
    threads = np.array([1, 2, 4, 8, 16])

#     boost = np.array([[269.324,124.678,67.460,35.700,18.003],
# [131.495,65.647,36.782,20.286,10.890],
# [66.814,36.521,21.548,11.698,7.158],
# [39.153,26.375,16.598,7.974,5.779],
# [34.167,26.067,17.014,6.613,5.764]])

    boost = np.array([[370.988,216.786,147.470,121.557,120.679],
[185.553,118.324,79.342,71.343,69.341],
[100.560,64.415,45.323,36.983,35.013],
[65.346,40.422,26.983,23.353,20.940],
[34.452,23.124,14.643,13.012,11.159]])
    boost = boost.T


#     boost = np.array([[48.056,23.525,12.852,8.055,5.644],
# [24.182,12.943,9.080,6.593,4.823],
# [13.626,8.131,5.620,4.639,4.058],
# [8.331,5.813,5.053,4.133,3.653],
# [6.150,5.117,4.456,3.678,2.968]])
#     boost = np.array([[281.437,137.551,61.090,30.058,16.856],
# [159.039,66.207,33.952,19.772,11.913],
# [90.105,46.236,21.548,13.165,8.532],
# [62.896,37.006,13.807,10.348,7.359],
# [53.235,28.961,12.879,8.652,6.123]])
    serial = boost[0][0]
    boost = np.array([[serial/boost[i][j] for j in range(len(boost[i]))] for i in range(len(boost))])
    # fault = (np.loadtxt("faults_implicit.txt"))
    
    fig = plt.figure("Boost")

    ax = fig.add_subplot(111, projection = '3d')

    x, y = np.meshgrid(nodes, threads)
    print (x)
    print (y)
    print(boost)


    ax.plot_surface(x, y, boost, alpha=0.9, rstride=1, cstride=1, linewidth=0.5, cmap='jet')
    # ax.scatter(x, y, fault)

    ax.set_xlabel('количество узлов')
    ax.set_ylabel('количество потоков')
    ax.set_zlabel('ускорение')
    ax.set_zbound(0, 50)

    # ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: '1e-{:d}'.format(int(x))))

    ax.view_init(azim=-120)
    plt.savefig('plot.png', transparent=True)
    plt.show()

if __name__ == '__main__':
    main()