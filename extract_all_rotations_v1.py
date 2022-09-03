from TrajectoryAnalyzerQuat import *
import sys

DIR = '../../../kjun_temporary_aimds/19_02_density_18/800K'
path = [x for x in os.listdir(DIR) if x.startswith('run')]
path.sort()
full_paths = ['{}/{}' .format(DIR, x) for x in path]

for i in full_paths:
    print(i)
result_rot = RotationAnalyzer.from_paths(full_paths, 'P', 800, step_skip=10, n_process=16)
result_rot.count_rot_from_graph_each_time_atomwise(n_process=16)#, radian)
result_rot.count_rot_from_graph_from_init_atomwise(n_process=16)#, radian)
result_rot.export_rotation_analysis(DIR)