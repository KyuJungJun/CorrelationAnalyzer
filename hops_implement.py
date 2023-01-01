from TrajectoryAnalyzerQuat import *
import sys
import argparse
from multiprocessing import set_start_method


DIR = '/Users/kjun/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Qualtimehome/Volumes/qualtime/Lacie_HDD_5TB/Blee/blee_backup/LPS_aimd/kjun_temporary_aimds/00_crystal_beta_314_try/hop_debug/1000K'
path = [x for x in os.listdir(DIR) if x.startswith('run')]
path.sort()
full_paths = ['{}/{}' .format(DIR, x) for x in path]

result_rot = HopAnalyzer.from_paths(full_paths[:4], 'Li', 1000, step_skip=10, n_process=16)
result_rot.count_hop_from_graph_each_time_atomwise(n_process=16, split=7500, part_index=int(0))#, radian)
result_rot.plot_hops(mode='hop_each_time')