# -*- coding: utf-8 -*-

import argparse as ap
import numpy as np

from matplotlib import pyplot, patches
from water_radius import parseTrajectories, parseResidues
from water_tracking import trackWaters


def parseArgs():
    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required=True, metavar="FILE", type=str, nargs='*', help="path to trajectory files")
    required.add_argument("-w", "--waters", metavar="CHAIN:ID", type=str, nargs='*', help="list of water ids", default=[])
    parser._action_groups.append(optional)
    args = parser.parse_args()

    trajectories = parseTrajectories(args.input, parser)
    waters = parseResidues(args.waters)

    return trajectories, waters


def main():
    trajectories, waters = parseArgs()
    num_waters = len(waters)
    print "Water shifts of {} molecule{} are going to be analyzed".format(num_waters, ["","s"][num_waters > 1])

    water_shifts = {}
    for water in waters:
        water_shifts[water[0] + water[1]] = []

    for trajectory in trajectories:
        water_tracking = trackWaters([trajectory, ], waters)
        for water, positions in water_tracking.iteritems():
            for index, position in enumerate(positions[:-1]):
                dist = np.linalg.norm([float(i) - j for i, j in zip(position, positions[index+1])])
                water_shifts[water].append(dist)

    for water, shifts in water_shifts.iteritems():
        print "Water {}:".format(water)
        print " - Mean shift: {} A".format(np.mean(shifts))
        print " - Variance:   {} A".format(np.var(shifts))
            

if __name__ == "__main__":
    main()
