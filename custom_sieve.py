#!/usr/bin/python3

import pandas as pd
import argparse as ap
import os
import glob


ACCEPTED_STEPS_COL = 'numberOfAcceptedPeleSteps'
WATER_DISTANCE_COL = 'COM DISTANCE'
BINDING_ENERGY_COL = 'Binding Energy'
CURRENT_ENERGY_COL = 'currentEnergy'
TRAJECTORY_NUM_COL = 'Trajectory Number'
INITIAL_STRUCT_COL = 'Initial Structure'
RMSD_DEVIATION_COL = 'proteinLigandDistance'
MAXIMUM_ACCEPTED_WATER_DISTANCE = 4.0
MIN_ACCEPTED_RMSD = 1.0


def parseArgs():
    working_dir = os.getcwd()

    parser = ap.ArgumentParser()
    parser.add_argument("-i", metavar="PATH", type=str, help="Path to PELE output files", default=working_dir)
    parser.add_argument("-o", metavar="PATH", type=str, help="Output path", default=working_dir)
    parser.add_argument("-d", metavar="FLOAT", type=float, help="Maximum accepted water distance", default=MAXIMUM_ACCEPTED_WATER_DISTANCE)
    parser.add_argument("-s", metavar="INT", type=int, help="Filter the results by an initial structure", default=None)
    args = parser.parse_args()

    in_path =  os.path.abspath(args.i)
    out_path =  os.path.abspath(args.o)
    return in_path, out_path, args.d, args.s


def getAllPeleReports(path):
    reports = glob.glob(os.path.join(path, "*report*"))
    if (len(reports) == 0):
        raise NameError('No Pele reports found in the input path') 
    return reports


def linkTrajectoriesWithSameStartingPoint(parsed_reports):
    trajectories = {}
    for index, row in parsed_reports.iterrows():
        if row[ACCEPTED_STEPS_COL] != 0:
            parsed_reports.loc[index, INITIAL_STRUCT_COL] = trajectories[initial_energy]
            continue
        initial_energy = row[CURRENT_ENERGY_COL]
        if initial_energy not in trajectories:
            trajectories[initial_energy] = len(trajectories) + 1
        parsed_reports.loc[index, INITIAL_STRUCT_COL] = trajectories[initial_energy]
    return parsed_reports


def parsePeleReport(path, report_id):
    parsed_report = pd.read_csv(path, sep='    ', engine='python')
    parsed_report[TRAJECTORY_NUM_COL] = report_id
    return parsed_report


def getWaterMediatedStructures(parsed_reports, accepted_wat_dist):
    all_reports = pd.concat(parsed_reports)
    all_trajectories = all_reports.loc[:, [TRAJECTORY_NUM_COL,
                                           ACCEPTED_STEPS_COL,
                                           WATER_DISTANCE_COL,
                                           RMSD_DEVIATION_COL,
                                           BINDING_ENERGY_COL,
                                           CURRENT_ENERGY_COL]]
    linked_trajectories = linkTrajectoriesWithSameStartingPoint(all_trajectories)
    selected_steps = linked_trajectories[linked_trajectories[WATER_DISTANCE_COL] < accepted_wat_dist]    

    return selected_steps


def sortReportsBy(parsed_reports, column, criteria='min'):
    if criteria is 'min':
        results = parsed_reports.sort_values(by=[column])
    elif criteria is 'max':
        results = parsed_reports.sort_values(by=[column], ascending=True)
    else:
        raise AttributeError('Unknown criteria')

    return results


def filterReportsBy(parsed_reports, column, value):
    results = parsed_reports[parsed_reports[column] == value]
    return results


def filterReportsByRMSD(parsed_reports):
    results = parsed_reports[parsed_reports[RMSD_DEVIATION_COL] <= MIN_ACCEPTED_RMSD]
    return results


def renameDataFrame(parsed_reports):
    results = parsed_reports.reset_index(drop=True)
    results = results.rename(index=str, columns={'Trajectory Number': 'Traj.',
                                                 'numberOfAcceptedPeleSteps': 'Step',
                                                 'COM DISTANCE': 'Wat. dist.',
                                                 'proteinLigandDistance': 'RMSD',
                                                 'Binding Energy': 'Bin. E.',
                                                 'currentEnergy': 'Tot. E.',
                                                 'Initial Structure': 'Ini. Struct.'})
    return results


def main(in_path, out_path, accepted_wat_dist, initial_struct):
    pele_reports = getAllPeleReports(in_path)

    parsed_reports = {}
    for report in pele_reports:
        report_id = os.path.basename(report).split("_")[-1]
        parsed_reports[report_id] = parsePeleReport(report, report_id)

    best_reports = getWaterMediatedStructures(parsed_reports, accepted_wat_dist)
    sorted_reports = sortReportsBy(best_reports, BINDING_ENERGY_COL)
    filtered_by_rms_reports = filterReportsByRMSD(sorted_reports)

    if initial_struct is not None:
        filtered_reports = filterReportsBy(filtered_by_rms_reports, INITIAL_STRUCT_COL, initial_struct)

    else:
        filtered_reports = filtered_by_rms_reports

    final_data_frame = renameDataFrame(filtered_reports)

    print(final_data_frame)
    

if __name__ == "__main__":
    in_path, out_path, accepted_wat_dist, initial_struct = parseArgs()
    main(in_path, out_path, accepted_wat_dist, initial_struct)
