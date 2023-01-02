from TrajectoryAnalyzerQuat import *
import sys
import argparse
from multiprocessing import set_start_method

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument("part_index")
    parser.add_argument("temperature")
    args = parser.parse_args()
    #print("Analyzing part : {}".format(args.part_index))
    print("Temperature: {}".format(args.temperature))
    DIR = '.'
    try:
        set_start_method('forkserver')
    except RuntimeError:
        pass
    DIR = '.'
    path = [x for x in os.listdir(DIR) if x.startswith('run')]
    path.sort()
    full_paths = ['{}/{}' .format(DIR, x) for x in path]
    result_rot = RotationAnalyzer.from_base(DIR, 'P', int(args.temperature), step_skip=10, n_process=16)
    result_rot.export_rot_graph_only()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=7500, part_index=0)#, radian)
    result_rot.export_rotation_analysis()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=7500, part_index=1)#, radian)
    result_rot.export_rotation_analysis()
    result_rot = RotationAnalyzer.from_many_npys(DIR)
    result_rot.plot_rotations('each_time')
    result_rot.count_rotations_single_dt(1000, 1) # 2 : hop distance
    result_rot.count_rotations_single_dt(1000, 20) # 2 : hop distance