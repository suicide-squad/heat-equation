import numpy as np
import matplotlib.pyplot as plt

import sys
import os

def main():
	if len(sys.argv) != 2:
		print("Не верное количество параметров! Укажите расположение "
              "файла относительно папки result")
		return

	pathSetting = os.path.join(os.pardir, "initial", "setting.ini")

	with open(pathSetting) as file:
		setting = {line.split('=')[0]: float(line.split('=')[1]) for line in file}

	xStart = setting['XSTART']
	xEnd = setting['XEND']
	NX = setting['NX']

	x = np.linspace(xStart, xEnd, NX, dtype = float)

	fileStart = os.path.join(os.pardir, "initial", "function.txt")
	fileFinish = os.path.join(os.pardir, "result", sys.argv[1])

	try:
		yStart = np.loadtxt(fileStart)[:len(x)]
		yFinish = np.loadtxt(fileFinish)

		assert len(x) == len(yFinish)

		# Рисование графиков
		plt.plot(x, yStart, label ='start time')
		plt.plot(x, yFinish, label ='end time')

		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
				   ncol=2, mode="expand", borderaxespad=0.)
		plt.xlabel('x', fontsize=14)
		plt.ylabel('U(x)', fontsize=14)
		plt.grid(True)

		plt.show()

	except AssertionError:
		print('НЕ совпадают размерности!')
	except FileNotFoundError:
		print('НЕ найден файл!')

if __name__ == '__main__':
	main()
