from sets import Set
from matplotlib import colors, cm

import chimera
from Shape.shapecmd import sphere_shape


COORDINATES_FILE = ".coordinates_file.tmp"


def plotPosition(coordinates_file):
	with open(coordinates_file, 'r') as cf:
		names = {}
		points = []
		reference = cf.readline().strip()
		for line in cf:
			name = line.split()[0]
			coordinates = line.split()[1:]
			if name not in names:
				names[name] = len(names)
			points.append((coordinates, names[name], name))

		color_map = cm.jet
		categories = len(names)

		for point in points:
			sphere_shape(radius=0.5,
						 divisions=10,
						 center=", ".join(point[0]),
						 color=[i * j for i, j in zip(color_map(point[1] * (color_map.N - 1) / (categories - 1)), [1, 1, 1, 0.3])],
						 modelId=point[1],
						 modelName=point[2])

		print "open %s" % reference
		chimera.runCommand("open {}".format(reference))


def main():
	plotPosition(COORDINATES_FILE)


if __name__ == "__main__":
	main()