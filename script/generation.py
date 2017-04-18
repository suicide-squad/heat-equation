
import numpy as np
from math import cos, pi
import os

def function(x, y, z):
    return cos(x * pi) if -0.5 < x < 0.5 else 0

def main():
    pathSetting = os.path.join(os.pardir, "initial", "setting.ini")
    with open(pathSetting) as file:
        setting = {line.split('=')[0] : float(line.split('=')[1]) for line in file}

    xStart = setting['XSTART']
    xEnd = setting['XEND']
    NX = setting['NX']

    yStart = setting['YSTART']
    yEnd = setting['YEND']
    NY = setting['NY']

    zStart = setting['ZSTART']
    zEnd = setting['ZEND']
    NZ = setting['NZ']

    X = np.linspace(xStart, xEnd, NX, dtype = float)
    Y = np.linspace(yStart, yEnd, NY, dtype = float)
    Z = np.linspace(zStart, zEnd, NZ, dtype = float)

    U = np.array([function(x, y, z) for z in Z for y in Y for x in X])

    pathOutput = os.path.join(os.pardir, "initial", 'function.txt')
    np.savetxt(pathOutput, U, fmt='%.15e',delimiter='\n')

if __name__ == "__main__":
    print("Start")
    main()
print("Finish")