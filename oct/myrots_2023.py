from TrajectoryAnalyzerOct import *
import sys
import argparse
from multiprocessing import set_start_method
import glob2

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

    DIR = '..'
    with open("trajectory.pickle", "wb") as f:
        structures = pickle.load(f)
    base_path = os.getcwd()+"/.."
    result_rot = RotationAnalyzer.from_structures(structures, base_path, 'Ta', ['O', 'Cl'], int(args.temperature), step_skip=10, n_process=16)
    result_rot.export_rot_graph_only()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=5000, part_index=0)#, radian)
    result_rot.export_rotation_analysis()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=5000, part_index=1)#, radian)
    result_rot.export_rotation_analysis()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=5000, part_index=3)#, radian)
    result_rot.export_rotation_analysis()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=5000, part_index=4)#, radian)
    result_rot.export_rotation_analysis()
    result_rot.count_rot_from_graph_each_time_atomwise(n_process=16, split=5000, part_index=5)#, radian)
    result_rot.export_rotation_analysis()

    result_rot = RotationAnalyzer.from_many_npys(DIR)
    print(result_rot.base)
    print(result_rot.rotations_each_time)
    result_rot.plot_rotations('each_time')