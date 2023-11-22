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
    result_hop = HopAnalyzer.from_base(DIR, 'Li', int(args.temperature), step_skip=10, n_process=16)
    result_hop.export_hop_graph_only()
    result_hop.count_hop_from_graph_each_time_atomwise(n_process=16, split=7500, part_index= 0) #int(args.part_index))#, radian)
    result_hop.export_hop_analysis()
    result_hop.count_hop_from_graph_each_time_atomwise(n_process=16, split=7500, part_index= 1)
    result_hop.export_hop_analysis()
    result_hop.count_hop_from_graph_each_time_atomwise(n_process=16, split=7500, part_index= 2)
    result_hop.export_hop_analysis()
    result_hop = HopAnalyzer.from_many_npys(DIR)
    result_hop.plot_hops('hop_each_time')
    result_hop.count_hops_single_dt(1000, 3) # 2 : hop distance
    result_hop.count_hops_single_dt(1000, 1) # 2 : hop distance