import argparse as ap
import os
import numpy as np
from matplotlib import pyplot


def parseArgs():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, metavar="FILE", type=str, help="path to input file")
    args = parser.parse_args()

    input_file =  os.path.abspath(args.input)
    if not os.path.exists(input_file):
        print "Error: path to reference \'", input_file, "\' not found."
        parser.print_help()
        exit(1)

    return input_file


def getData(coordinates_file):
	with open(coordinates_file, 'r') as cf:
		data = {}
		reference = cf.readline().strip()
		for line in cf:
			name = line.split()[0]
			coordinates = line.split()[1:]
			if name not in data:
				data[name] = []
			data[name].append(coordinates)

	return data


def analyzeData(data):
	for water in data:
		print("\nWater {}".format(water))
		x = [float(i[0]) for i in data[water]]
		y = [float(i[1]) for i in data[water]]
		z = [float(i[2]) for i in data[water]]
		x_mean = np.mean(x)
		y_mean = np.mean(y)
		z_mean = np.mean(z)
		x_var = np.var(x)
		y_var = np.var(y)
		z_var = np.var(z)

		dist = []
		for ind, point in enumerate(data[water]):
			current_dist = np.linalg.norm([float(i) - j for i, j in zip(point, (x_mean, y_mean, z_mean))])
			dist.append(current_dist)
			data[water][ind] = [x[ind], y[ind], z[ind], current_dist]
		central_dist = np.mean(dist)

		print("x:\n\t- mean     = {: 7.3f}\n\t- variance = {: 7.3f}".format(x_mean, x_var))
		print("y:\n\t- mean     = {: 7.3f}\n\t- variance = {: 7.3f}".format(y_mean, y_var))
		print("z:\n\t- mean     = {: 7.3f}\n\t- variance = {: 7.3f}".format(z_mean, z_var))
		print("central point:\n\t({: 7.3f},{: 7.3f},{: 7.3f})".format(x_mean, y_mean, z_mean))
		print("average distance from central point:\n\t{: 7.3f}".format(central_dist))

	return data


def plotData(data):
	pyplot.boxplot([[i[3] for i in data[j]] for j in data], labels=[i for i in data])
	pyplot.show()


def main():
	input_file = parseArgs()
	data = getData(input_file)
	data = analyzeData(data)
	plotData(data)


if __name__ == "__main__":
	main()