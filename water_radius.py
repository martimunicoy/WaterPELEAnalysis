# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import argparse as ap
import os
import glob
import numpy as np
from sets import Set
import copy
from matplotlib import pyplot, rc
from matplotlib import patches
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
        print "Warning: list of water ids is empty. No correct water ids were detected."

    return waters


def parseTrajectories(trajectories_to_parse, parser):
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
    required.add_argument("-r", "--ref", required=True, metavar="FILE", type=str, help="path to reference structure file")
    required.add_argument("-w", "--waters", metavar="CHAIN:ID", type=str, nargs='*', help="list of water ids", default=[])
    required.add_argument("-i", "--input", required=True, metavar="FILE", type=str, nargs='*', help="path to trajectory files")
    optional.add_argument("-R", "--radius", metavar="FLOAT", type=float, help="radius of the sphere to look for waters", default=1.5)
    optional.add_argument("-X", "--xaxis", metavar="INTEGER [METRIC]", type=str, nargs='*', help="column number and metric to plot on the X axis", default=None)
    optional.add_argument("-Y", "--yaxis", metavar="INTEGER [METRIC]", type=str, nargs='*', help="column number and metric to plot on the Y axis", default=None)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    reference =  os.path.abspath(args.ref)
    if not os.path.exists(reference):
        print "Error: path to reference \'", reference, "\' not found."
        parser.print_help()
        exit(1)

    waters = parseResidues(args.waters)

    trajectories = parseTrajectories(args.input, parser)

    radius = args.radius

    x_data = args.xaxis
    y_data = args.yaxis

    return reference, waters, trajectories, radius, x_data, y_data


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
    matchs = {}

    for trajectory in trajectories:
        traj_directory = os.path.dirname(trajectory)
        traj_number = os.path.basename(trajectory).split('_')[-1].split('.')[0]

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

                coordinates = fields[6:9]
                results[model][fields[4] + fields[5]] = waterInSphere(coordinates, water_locations, radius)

        matchs[traj_directory, traj_number] = {}

        for i in xrange(num_waters + 1):
            matchs[traj_directory, traj_number][i] = []

        for model, water_matchs in results.iteritems():
            sorted_waters = sorted(water_matchs, key=lambda k: len(water_matchs[k]), reverse=True)
            if len(water_matchs[sorted_waters[0]]) == 0:
                continue

            match_set = []
            for water in sorted_waters:
                for match in water_matchs[water]:
                    if match not in match_set:
                        match_set.append(match)
                        break

            matchs[traj_directory, traj_number][len(match_set)].append(model)

    return matchs


def parseAxisData(axis_data):
    if axis_data is None:
        return ([None, ] , None)
    else:
        try:
            rows = [int(axis_data[0]), ]
        except ValueError:
            print("Warning: axis data not recognized: {}".format(axis_data))
            return ([None, ], None)
        if len(axis_data) == 1:
            return (rows, '?')
        elif len(axis_data) > 1:
            label_index = 1
            while axis_data[label_index] == "+":
                label_index += 2
                try:
                    rows.append(int(axis_data[2]))
                except (ValueError, IndexError):
                    print("Warning: axis data not recognized: {}".format(axis_data))
                    return ([None, ], None)
                if len(axis_data) == label_index:
                    return (rows, '?')
            if "energy" in axis_data[label_index].lower():
                label = axis_data[label_index] + " ($kcal/mol$)"
            elif "energies" in axis_data[label_index].lower():
                label = axis_data[label_index] + " ($kcal/mol$)"
            elif "distance" in axis_data[label_index].lower():
                label = axis_data[label_index] + " ($\AA$)"
            elif "rmsd" in axis_data[label_index].lower():
                label = axis_data[label_index] + " ($\AA$)"
            else:
                label = axis_data[label_index]
            return (rows, label)

        print("Warning: axis data not recognized: {}".format(axis_data))
        return ([None, ], None)


def scatterPlot(matchs, x_rows=[None, ], y_rows=[None, ], x_name=None, y_name=None):
    x_values = []
    y_values = []
    labels = []
    annotations = []

    if None in x_rows:
        x_rows = [7, ]
        x_name = "RMSD ($\AA$)"
    if None in y_rows:
        y_rows = [5, ]
        y_name = "Energy ($kcal/mol$)"
    if x_name is None:
        x_name = '?'
    if y_name is None:
        y_name = '?'

    for traj_info, categories in matchs.iteritems():
        traj_directory, traj_number = traj_info
        labels_size = len(labels)

        report = traj_directory + "/report_" + traj_number

        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                total = 0.
                for x_row in x_rows:
                    total += float(line.split()[x_row - 1])
                x_values.append(total)
                total = 0.
                for y_row in y_rows:
                    total += float(line.split()[y_row - 1])
                y_values.append(total)
                epoch = traj_directory.split('/')[-1]
                if not epoch.isdigit():
                    epoch = '0'
                    
                annotations.append("Epoch: " + epoch + "\n" + "Trajectory: " + traj_number + "\n" + "Model: " + str(i + 1))

                labels.append(0)

        for category, models in categories.iteritems():
            for model in models:
                labels[labels_size + model - 1] = category + 1

    norm = pyplot.Normalize(1, max(matchs.values()[0]) + 1)
    cmap = pyplot.cm.RdYlGn

    fig, ax = pyplot.subplots()
    sc = pyplot.scatter(x_values, y_values, c=labels, cmap=cmap, norm=norm, alpha=0.6)

    ax.margins(0.05)
    ax.set_facecolor('lightgray')
    pyplot.ylabel(y_name)
    pyplot.xlabel(x_name)

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    patches_list = [patches.Patch(color=cmap(norm(1)), label='No matches', alpha=0.6), ]
    for i in xrange(max(matchs.values()[0])):
        if i == 0:
            match_str = "match"
        else:
            match_str = "matches"
        patches_list.append(patches.Patch(color=cmap(norm(i + 2)), label='{} {}'.format(i + 1, match_str), alpha=0.6))

    ax.legend(handles=patches_list)

    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        annot.set_text(annotations[int(ind["ind"][0])])
        annot.get_bbox_patch().set_facecolor(cmap(norm(labels[ind["ind"][0]])))

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    pyplot.show()

 
def main():
    reference, waters, trajectories, radius, x_data, y_data = parseArgs()
    num_waters = len(waters)
    water_locations = getWaterReferenceLocations(reference, waters)
    matchs = findWaterMatchs(trajectories, waters, water_locations, radius, num_waters)

    x_rows, x_name = parseAxisData(x_data)
    y_rows, y_name = parseAxisData(y_data)

    scatterPlot(matchs, x_rows=x_rows, y_rows=y_rows, x_name=x_name, y_name=y_name)


if __name__ == "__main__":
    main()
