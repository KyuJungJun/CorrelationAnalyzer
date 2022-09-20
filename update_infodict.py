from TrajectoryAnalyzerQuat import *
import sys
import argparse
from multiprocessing import set_start_method
DIR = '.'
rot = RotationAnalyzer.from_many_npys(DIR, read_structures=False)
if 'init_structure' not in rot.info_dict.keys():
	print("Extracting rot graph")
	path = [x for x in os.listdir(DIR) if x.startswith('run')]
	path.sort()
	full_paths = ['{}/{}' .format(DIR, x) for x in path]
	s = Vasprun(full_paths[0]+'/vasprun.xml').structures[0]
	print(s)
	rot.info_dict['init_structure'] = s
	rot.export_rot_graph_only(DIR)
	print('exported again')
rot = RotationAnalyzer.from_many_npys(DIR, read_structures=False)
print(rot.info_dict['init_structure'])
print(rot.rotations_each_time.shape)