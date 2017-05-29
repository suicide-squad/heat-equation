import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import matplotlib.pyplot as plt
import numpy as np
import os

def buildPlot3D(ax, title, u, rstride, cstride, usmin, usmax,
                 x, y, z, xmax, xmin, ymax, ymin, zmax, zmin):
    norm = colors.Normalize(vmin=usmin, vmax=usmax, clip=True)
    alpha = 0.8
    linewidth = 0.5
    antialiased = True
    shade = False

    X, Y = np.meshgrid(x, y)

    ax.plot_surface(X, Y, zmax, alpha=alpha, facecolors=cm.jet(norm(u[:][:][-1])),
                    rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)
    ax.plot_surface(X, Y, zmin, alpha=alpha, facecolors=cm.jet(norm(u[:][:][0])),
                    rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)

    X, Z = np.meshgrid(x, z)

    ax.plot_surface(X, ymax, Z, alpha=alpha, facecolors=cm.jet(norm(u[:][-1][:])),
                    rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)
    ax.plot_surface(X, ymin, Z, alpha=0.7, facecolors=cm.jet(norm(u[:][0][:])),
                    rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)

    Y, Z = np.meshgrid(y, z)

    ax.plot_surface(xmax, Y, Z, alpha=alpha, facecolors=cm.jet(norm(u[-1][:][:])),
                    rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)
    ax.plot_surface(xmin, Y, Z, alpha=alpha, facecolors=cm.jet(norm(u[0][:][:])),
                               rstride=rstride, cstride=cstride, linewidth=linewidth,
                    antialiased=antialiased, shade=shade)

    # plt.title(title)

    # ax.set_xlim(xmin + 1, xmax - 1)
    # ax.set_ylim(ymin + 1, ymax - 1)
    # ax.set_zlim(zmin + 1, zmax - 1)


    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

def main():
    if len(sys.argv) != 4:
        print("Не верное количество параметров! Укажите расположение "
              "файла настроек, начальную функцию относительно папки initial"
              " и файл результирующей функции относительно папки result.")
        return

    pathSetting = os.path.join(os.pardir, "initial", sys.argv[1])
    with open(pathSetting) as file:
        setting = {line.split('=')[0] : float(line.split('=')[1]) for line in file}

    xmin = setting['XSTART']
    xmax = setting['XEND']
    NX = int(setting['NX'])

    ymin = setting['YSTART']
    ymax = setting['YEND']
    NY = int(setting['NY'])

    zmin = setting['ZSTART']
    zmax = setting['ZEND']
    NZ = int(setting['NZ'])

    tstart = 'START {:.0e} sec'.format(setting['TSTART'])
    tfinish = 'FINISH {:.0e} sec'.format(setting['TFINISH'])

    x = np.linspace(xmin, xmax, NX, dtype = float)
    y = np.linspace(ymin, ymax, NY, dtype = float)
    z = np.linspace(zmin, zmax, NZ, dtype = float)

    pathStart = os.path.join(os.pardir, "initial", sys.argv[2])
    us = np.loadtxt(pathStart)
    us.shape = (NZ, NY, NX)

    usmin = -1
    usmax = 1

    cstride = max(NX // 20, 1)
    rstride = max(NY // 20, 1)
    # cstride = 1
    # rstride = 1
    print(rstride, cstride)

    fig = plt.figure("HEAD EQUATION")

    ax1 = fig.add_subplot(121, projection = '3d')
    buildPlot3D(ax1, tstart, us, rstride, cstride, usmin, usmax, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax)

    pathFinish = os.path.join(os.pardir, "result", sys.argv[3])
    uf = np.loadtxt(pathFinish)
    uf.shape = (NZ, NY, NX)

    ax2 = fig.add_subplot(122, projection = '3d')
    buildPlot3D(ax2, tfinish, uf, rstride, cstride, usmin, usmax, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax)

    m = cm.ScalarMappable(cmap=cm.jet)
    head = np.linspace(-1.0, 1.0, 100)
    m.set_array(head)

    cax = fig.add_axes([0.12, 0.05, 0.78, 0.02])
    fig.colorbar(m, cax=cax, orientation='horizontal')

    plt.savefig('plot3D.png', transparent=True)
    plt.show()

if __name__ == '__main__':
    main()