from TrajectoryAnalyzerQuat import *
import sys
import argparse
from multiprocessing import set_start_method
DIR = '.'
rot = RotationAnalyzer.from_many_npys(DIR, read_structures=False)
ignored_angles = 15
min_dt = 100
rot_list = rot.count_rotations_from_each_time(ignored_angles=ignored_angles,min_dt=min_dt)
for index, i  in enumerate(rot_list):
	print(index)
	for rot_item in i:
		print(rot_item.time, rot_item.angle, rot_item.duration)
import pickle
with open('counter_{}_{}.pkl'.format(ignored_angles, min_dt), 'wb') as f:
	pickle.dump(rot_list, f)