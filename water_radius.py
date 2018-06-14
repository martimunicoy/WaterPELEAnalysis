# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import argparse as ap
import os
import glob
import numpy as np
from sets import Set
import copy
from matplotlib import pyplot, rc
import pandas as pd


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


def parseArgs():
    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-r", "--ref", required=True, metavar="PATH", type=str, help="path to reference structure")
    required.add_argument("-w", "--waters", required=True, metavar="CHAIN:ID", type=str, nargs='*', help="list of water ids")
    required.add_argument("-i", "--input", required=True, metavar="PATH", type=str, nargs='*', help="path to trajectory files")
    optional.add_argument("-R", "--radius", metavar="FLOAT", type=float, help="radius of the sphere to look for waters", default=1.5)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reference =  os.path.abspath(args.ref)
    if not os.path.exists(reference):
        print "Error: path to reference \'", reference, "\' not found."
        parser.print_help()
        exit(1)

    waters = parseResidues(args.waters)

    trajectories = parseTrajectories(args.input)

    radius = args.radius

    return reference, waters, trajectories, radius


def getWaterReferenceLocations(reference, waters):
    water_locations = []
    waters_list =  copy.copy(waters)
    with open(reference, "r") as ref_pdb:
        for line in ref_pdb:
            if not line.startswith('HETATM'):
                continue
            fields = line.split()
            if fields[3] != 'HOH' and fields[2] != 'O':
                continue
            for i, water in enumerate(waters_list):
                chain, residue_id = water
                if fields[4] == chain and fields[5] == residue_id:
                    water_locations.append([j for j in fields[6:9]])
                    del(waters_list[i])
                    break

    if len(waters_list) != 0:
        print "Warning: the following water residues could not be found in the reference structure:"
        for water in waters_list:
            print water

    return water_locations


def waterInSphere(coordinates, water_locations, radius):
    squared_radius = pow(radius, 2)
    matchs = []
    for k, water_location in enumerate(water_locations):
        squared_distance = sum([pow(float(i) - float(j), 2) for i, j in zip(coordinates, water_location)])
        if squared_distance < squared_radius:
            matchs.append(k)
    return matchs


def findWaterMatchs(trajectories, waters, water_locations, radius, num_waters):
    bests = []
    for trajectory in trajectories:
        with open(trajectory, "r") as pdb_file:
            results = {}
            model = int(pdb_file.readline().split()[1])
            results[model] = {}
            for line in pdb_file:
                if line.startswith("MODEL"):
                    model += 1
                    results[model] = {}
                    continue
                if not line.startswith('HETATM'):
                    continue
                fields = line.split()
                if fields[3] != 'HOH' or fields[2] != 'OW':
                    continue
                if (fields[4], fields[5]) in waters:
                    coordinates = fields[6:9]
                    results[model][fields[4] + fields[5]] = waterInSphere(coordinates, water_locations, radius)

        if len(results[model]) != num_waters:
            print "Error: some waters are missing in trajectory files:"
            for water in waters:
                if water[0] + water[1] not in results[model]:
                    print "Water: ", water[0], ":", water[1]
            exit(1)

        for model, water_matchs in results.iteritems():
            sorted_waters = sorted(water_matchs, key=lambda k: len(water_matchs[k]))
            if len(water_matchs[sorted_waters[0]]) == 0:
                continue

            match_set = []
            for water in sorted_waters:
                for match in water_matchs[water]:
                    if match not in match_set:
                        match_set.append(match)
                        break

            if len(match_set) == num_waters:
                bests.append([trajectory, model])

    return bests


def addAttributesFromReport(bests):
    for best in bests:
        report_file = os.path.dirname(best[0]) + '/report_' + os.path.basename(best[0]).split('_')[-1].split('.')[0]
        with open(report_file, 'r') as report_file:
            for i, line in enumerate(report_file):
                if i == best[1]:
                    best.append(float(line.split()[3]))
                    best.append(float(line.split()[4]))
                    best.append(float(line.split()[6]))


def scatterPlot(trajectories, bests, x_row=6, y_row=4):
    reports = []
    for trajectory in trajectories:
        for report in glob.glob(os.path.dirname(trajectory) + '/report_' + os.path.basename(trajectory).split('_')[-1].split('.')[0]):
            reports.append(report)

    bests_ids = []
    for best in bests:
        bests_ids.append((os.path.dirname(best[0]), os.path.basename(best[0]).split('_')[-1].split('.')[0], best[1]))

    x_values = []
    y_values = []
    categories = []

    for report in reports:
        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                x_values.append(float(line.split()[x_row]))
                y_values.append(float(line.split()[y_row]))
                categories.append((os.path.dirname(report), os.path.basename(report).split('_')[-1].split('.')[0], i + 1) in bests_ids)

    data = pd.DataFrame(dict(x=x_values, y=y_values, label=categories))
    groups = data.groupby('label')

    fig, ax = pyplot.subplots()
    ax.margins(0.05)
    pyplot.ylabel("Energy ($kcal/mol$)")
    pyplot.xlabel("RMSD ($\AA$)")

    for name, group in groups:
        if name:
            label = "Match"
            color = "green"
        else:
            label = "No match"
            color = "red"
        ax.plot(group.x, group.y, ms=7, linestyle='', marker='o', markeredgewidth=0.0, alpha=0.6, label=label, color=color)

    ax.legend()

    pyplot.show()

 
def main():
    reference, waters, trajectories, radius = parseArgs()
    num_waters = len(waters)
    water_locations = getWaterReferenceLocations(reference, waters)
    bests = findWaterMatchs(trajectories, waters, water_locations, radius, num_waters)
    addAttributesFromReport(bests)

    if len(bests) == 0:
        print "No significant matches were found :("
    else:
        for best in bests:
            print "trajectory: ", best[0].split('_')[-1].split('.')[0], " pele step: ", best[1] - 1, " model: ", best[1], " total energy: ", best[2], " binding energy: ", best[3], " RMSD: ", best[4]

    scatterPlot(trajectories, bests, x_row=6, y_row=4)


if __name__ == "__main__":
    main()
