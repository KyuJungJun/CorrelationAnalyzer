from TrajectoryAnalyzerQuat import *
import sys
import argparse

if __name__ == '__main__':

    print("Extracting rot graph")
    DIR = '.'
    path = [x for x in os.listdir(DIR) if x.startswith('run')]
    path.sort()
    full_paths = ['{}/{}' .format(DIR, x) for x in path]
    for i in full_paths:
        print(i)
    result_rot = RotationAnalyzer.from_paths(full_paths, 'P', 700, step_skip=10, n_process=16)
    result_rot.export_rot_graph_only('.')