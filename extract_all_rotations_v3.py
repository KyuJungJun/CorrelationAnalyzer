from TrajectoryAnalyzerQuat import *
import sys
import argparse
from multiprocessing import set_start_method

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("part_index")
    args = parser.parse_args()
    print("Analyzing part : {}".format(args.part_index))
    DIR = '.'
    try:
        set_start_method('forkserver')
    except RuntimeError:
        pass
    '''
    path = [x for x in os.listdir(DIR) if x.startswith('run')]
    path.sort()
    full_paths = ['{}/{}' .format(DIR, x) for x in path]

    for i in full_paths:
        print(i)
    result_rot = RotationAnalyzer.from_paths(full_paths, 'P', 700, step_skip=10, n_process=16)
    '''
    result_rot = RotationAnalyzer.from_rot_graph(DIR)
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=7500, part_index=int(args.part_index))#, radian)
    result_rot.count_rot_from_graph_from_init_atomwise(n_process=16, split=7500, part_index=int(args.part_index))# radian)
    result_rot.export_rotation_analysis(DIR)