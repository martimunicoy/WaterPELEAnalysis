# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import argparse as ap
import os
import glob
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from subprocess import call

FILENAME = "WaterTracking"
CHIMERA_PATH = "/home/municoy/.local/UCSF-Chimera64-1.12/bin/chimera"


def parseTrajectories(trajectories_to_parse):
    trajectories = []

    for trajectory_list in trajectories_to_parse:
        trajectories_found = glob.glob(trajectory_list)
        if len(trajectories_found) == 0:
            print "Warning: trajectory path \'", trajectory_list, "\' not found."
        for trajectory in glob.glob(trajectory_list):
            trajectories.append(trajectory)

    if len(trajectories) == 0:
        print "Error: list of trajectories is empty."
        parser.print_help() 
        exit(1)

    return trajectories


def parseResidues(residues_to_parse):
    waters = []

    for water_list in residues_to_parse:
        for water in water_list.split(','):
            water.strip()
            water_identifiers = water.split(':')
            if len(water_identifiers) == 2:
                chain, residue_id = water_identifiers
                waters.append((chain, residue_id))

    if len(waters) == 0:
        print "Error: list of water ids is empty. No correct water ids were detected."
        parser.print_help()
        exit(1)

    return waters


def parseArgs():
    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="PATH", type=str, nargs='*', help="path to trajectory files")
    required.add_argument("-w", "--waters", required=True, metavar="CHAIN:ID", type=str, nargs='*', help="list of water ids")
    required.add_argument("-r", "--ref", required=True, metavar="PATH", type=str, help="path to reference structure")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reference =  os.path.abspath(args.ref)
    if not os.path.exists(reference):
        print "Error: path to reference \'", reference, "\' not found."
        parser.print_help()
        exit(1)
    trajectories = parseTrajectories(args.input)
    waters = parseResidues(args.waters)

    return reference, trajectories, waters


def trackWaters(trajectories, waters):
    results = {}
    
    for water in waters:
        results[water[0] + water[1]] = []

    for trajectory in trajectories:
        with open(trajectory, 'r') as trajectory_file:
            for line in trajectory_file:
                if not line.startswith('HETATM'):
                    continue
                fields = line.split()
                if fields[3] != 'HOH' or fields[2] != 'OW':
                    continue
                if fields[4] + fields[5] in results:
                    coordinates = tuple(map(float, fields[6:9]))
                    results[fields[4] + fields[5]].append(coordinates)

    return results


def plotWaterTracking(data):
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    for water, coordinates in data.iteritems():
        ax.scatter([i[0] for i in coordinates], [j[1] for j in coordinates], [k[2] for k in coordinates])

    ax.set_xlabel('X ($\AA$)')
    ax.set_ylabel('Y ($\AA$)')
    ax.set_zlabel('Z ($\AA$)')

    pyplot.show()


def saveTrackingToPDB(data, reference):
    filename_path = os.path.abspath(FILENAME + ".pdb")
    filename_dir = os.path.dirname(filename_path)
    filename_id = 0

    while os.path.exists(filename_path):
        filename_id += 1
        filename_path = filename_dir + '/' + FILENAME + "_" + str(filename_id) + ".pdb"

    with open(filename_path, 'w') as pdb_file:
        with open(reference, 'r') as ref_file:
            for line in ref_file:
                if (line.startswith("TITLE") or
                    line.startswith("MODEL") or
                    line.startswith("ATOM") or
                    line.startswith("TER") or
                    line.startswith("HETATM")):
                    pdb_file.write(line)
        index = 1
        pdb_file.write("TER\n")
        for water, coordinates in data.iteritems():
            for point in coordinates:
                pdb_file.write("HETATM ")
                pdb_file.write("%4i " %index)
                pdb_file.write(" WAT ")
                pdb_file.write("PNT ")
                pdb_file.write("%s " %water[0])
                pdb_file.write("%3d    " %float(water[1:]))
                pdb_file.write("% 7.3f " % (float(point[0])))
                pdb_file.write("% 7.3f " % (float(point[1])))
                pdb_file.write("% 7.3f " % (float(point[2])))
                pdb_file.write(" 1.00 ")
                pdb_file.write("31.67 ")
                pdb_file.write("           ")
                pdb_file.write("C  ")
                pdb_file.write("\n")
                index += 1
        pdb_file.write("END\n")
        pdb_file.write("\n")

    return filename_path


def saveCoordinatesFile(data, reference):
    filename_path = os.path.abspath(FILENAME + ".out")
    filename_dir = os.path.dirname(filename_path)
    filename_id = 0

    while os.path.exists(filename_path):
        filename_id += 1
        filename_path = filename_dir + '/' + FILENAME + "_" + str(filename_id) + ".out"

    with open(filename_path, 'w') as pdb_file:
        pdb_file.write(reference + '\n')
        for water, coordinates in data.iteritems():
            for point in coordinates:
                pdb_file.write("{} {: 7.3f} {: 7.3f} {: 7.3f}\n".format(water, float(point[0]), float(point[1]), float(point[2])))

    return filename_path

def main():
    reference, trajectories, waters = parseArgs()
    print "Tracking waters..."
    water_tracking = trackWaters(trajectories, waters)
    #plotWaterTracking(water_tracking)
    #filename_path = saveTrackingToPDB(water_tracking, reference)
    print "Saving coordinates..."
    filename_path = saveCoordinatesFile(water_tracking, reference)
    print "Coordinates saved at:", filename_path


if __name__ == "__main__":
    main()