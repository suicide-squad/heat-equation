import numpy as np
from math import cos, pi

def function(x, y, z):
    return cos(x * pi) if -0.5 < x < 0.5 else 0

with open('setting.ini') as file:
    setting = {line.split('=')[0] : float(line.split('=')[1]) for line in file}
    file.close()
print(setting)

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

U = np.array([function(x, y, z) for x in X for y in Y for z in Z])

np.savetxt('function.txt', U, delimiter='\n')
