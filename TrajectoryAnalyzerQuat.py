#!/usr/bin/env python
# coding: utf-8

# In[1]:
'''
CorrelationAnalyzer
Version: September 13nd, 2022 V0.3
Writer: KyuJung Jun (kjun@berkeley.edu), Dr. Byungju Lee (blee89@kist.re.kr)
'''

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
#from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder, find_rotation
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import PeriodicSite
from scipy.spatial.transform import Rotation as R
import os, sys, json, copy, time, itertools, multiprocessing, math
import numpy as np
from pymatgen.util.plotting import pretty_plot
import matplotlib.pyplot as plt
from pyquaternion import Quaternion
import scipy
import pickle
from multiprocessing import set_start_method
import sys

try:
    set_start_method('forkserver')
except RuntimeError:
    pass


# In[2]:
def angle_quat(x,y):
    print("Deprecated, do not use Quaternion.sym_distance")
    '''
    x, y: PyQuaternion objects
    TBD: check what is the difference in various distance metric for quaternions
    Note: This formulation is more numerically stable when performing iterative gradient descent
    on the Riemannian quaternion manifold.
    However, the distance between q and -q is equal to pi, rendering this formulation not useful
    for measuring rotation similarities when the samples are spread over a "solid" angle of more than pi/2 radians
    (the spread refers to quaternions as point samples on the unit hypersphere).
    '''
    return Quaternion.sym_distance(x,y)
"""
def abs_angle_quat(x,y):
    return Quaternion.distance(x,y)

"""
def abs_angle_quat(x,y):
    '''
    x, y: PyQuaternion objects
    returns the angle in radians
    Note: This function does not measure the distance on the hypersphere,
    but it takes into account the fact that q and -q encode the same rotation.
    It is thus a good indicator for rotation similarities.
    '''
    return Quaternion.absolute_distance(x,y)

def angle_two_vector(va, vb):
    cosine_angle = np.round(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)), decimals=5)
    return np.degrees(np.arccos(cosine_angle))

def produce_sorted_coords(neigbors_sets):
    unsorted_coords = neigbors_sets.neighb_coords - neigbors_sets.coords[0]
    index_store = []
    for atom in neigbors_sets.neighb_sites_and_indices:
        index_store.append(atom['index'])
    index_store.sort()
    sorted_index = []
    for ind in index_store:
        for i, atom in enumerate(neigbors_sets.neighb_sites_and_indices):
            if atom['index'] == ind:
                sorted_index.append(i)
    sorted_coords = []
    for i in sorted_index:
        sorted_coords.append(unsorted_coords[i])
    return np.array(sorted_coords)

def calc_quat(neigbors_sets1, neigbors_sets2):
    relative_coords1 = produce_sorted_coords(neigbors_sets1)
    relative_coords2 = produce_sorted_coords(neigbors_sets2)
    rot = find_rotation(relative_coords1, relative_coords2) # no permutation
    quat = R.from_matrix(rot).as_quat()
    Q = Quaternion(x=quat[0],y=quat[1],z=quat[2],w=quat[3])
    #angles = np.round(R.from_matrix(rot).as_euler('ZXZ', degrees=True), 1) 
    #angles = np.round(R.from_matrix(rot).as_quat())#('zxy', degrees=True), 1)
    return Q

def derive_neighbor_sets(strt, species, starting_index, num_indices, i, tot_num):
    print('Calculating {}-th structure out of {} structures' .format(i, tot_num))
    with HiddenPrints():
        LGF = LocalGeometryFinder()
        LGF.setup_structure(structure=strt)
        se = LGF.compute_structure_environments(only_atoms=[species], maximum_distance_factor=1.5, get_from_hints=True)
        ns_list = []
    for x in np.arange(starting_index, starting_index + num_indices):
        if 4 in se.neighbors_sets[x].keys():
            #if len(se.neighbors_sets[x][4])>1:
            #    print(se.neighbors_sets[x][4], len(se.neighbors_sets[x][4]))
            ns_list.append(se.neighbors_sets[x][4][0])
        else:
            ns_list.append(None)
     
    return ns_list

# In[3]:


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

'''
class DA_new(DiffusionAnalyzer):
    def __init__(self, structure, displacements, specie, temperature,
                 time_step, step_skip, smoothed="max", min_obs=30,
                 avg_nsteps=1000, lattices=None, disp_vs_time=None):
        super(DA_new, self).__init__(structure, displacements, specie, temperature,
                 time_step, step_skip, smoothed, min_obs,
                 avg_nsteps, lattices)
        self.disp_time = disp_vs_time

                                     
    @classmethod
    def from_structures(cls, structures, specie, temperature,
                        time_step, step_skip, initial_disp=None,
                        initial_structure=None, **kwargs):
        p, l = [], []
        for i, s in enumerate(structures):
            if i == 0:
                structure = s
            p.append(np.array(s.frac_coords)[:, None])
            l.append(s.lattice.matrix)
        if initial_structure is not None:
            p.insert(0, np.array(initial_structure.frac_coords)[:, None])
            l.insert(0, initial_structure.lattice.matrix)
        else:
            p.insert(0, p[0])
            l.insert(0, l[0])

        p = np.concatenate(p, axis=1)
        dp = p[:, 1:] - p[:, :-1]
        dp = dp - np.round(dp)
        dp_lattice = []
        for i in dp:
            dp_lattice.append([np.dot(d, m) for d, m in zip(i, l[1:])])
        f_disp = np.cumsum(dp, axis=1)
        c_disp = []
        for i in f_disp:
            c_disp.append([np.dot(d, m) for d, m in zip(i, l[1:])])
        disp = np.array(c_disp)

        # If is NVT-AIMD, clear lattice data.
        if np.array_equal(l[0], l[-1]):
            l = np.array([l[0]])
        else:
            l = np.array(l)
        if initial_disp is not None:
            disp += initial_disp[:, None, :]

        return cls(structure, disp, specie, temperature, time_step,
                   step_skip=step_skip, lattices=l, disp_vs_time=np.array(dp_lattice), **kwargs)

'''
# In[4]:


class Hop:
    def __init__(self, time, site, index, vector, hopping_length):
        self.time = time
        self.fin_site = site
        self.index = index
        self.vector = vector
        self.ini_site = copy.copy(self.fin_site)
        self.ini_site.coords = self.ini_site.coords - self.vector
        self.threshold_length = hopping_length
        
class HHCorrelation:
    def __init__(self, Hop1, Hop2, time_scale, length_scale):
        self.indice = [Hop1.index, Hop2.index]
        self.time_scale = time_scale
        self.length_scale = length_scale
        self.timestep = abs(Hop1.time - Hop2.time)
        
#         if Hop1.time <= Hop2.time:
#             self.distance = Hop1.fin_site.distance(Hop2.ini_site)
#         else:
#             self.distance = Hop1.ini_site.distance(Hop2.fin_site)
        self.distance = np.average([Hop1.ini_site.distance(Hop2.ini_site), Hop1.fin_site.distance(Hop2.fin_site)])
        self.angle = angle_two_vector(Hop1.vector, Hop2.vector)
        
              
class HoppingAnalyzer: 
    def __init__(self, structures, species, temperature, step_skip=10, time_step=2, hopping_length=3, text_input=None):
        self.species = Element(species)
        self.temperature = temperature
        self.step_skip = step_skip
        self.time_step = time_step
        self.hop_list = []
        
        if structures:
            self.dif = DA_new.from_structures(structures, species, temperature, self.time_step, self.step_skip, smoothed=False)
            self.disp_coord = []
            indices = []
            for i, k in enumerate(self.dif.disp_time):
                if structures[0].species[i] == Element(species):
                    self.disp_coord.append(k)
                    indices.append(i)
            self.disp_coord = np.array(self.disp_coord)
            starting_index = min(indices)
            print('The number of {} atoms: {}' .format(species, len(self.disp_coord)))

            for a, atom in enumerate(self.disp_coord):
                cumul_sum = np.zeros(3)
                for t, time in enumerate(atom):
                    cumul_sum = cumul_sum + time
                    if np.linalg.norm(cumul_sum) > hopping_length:
                        self.hop_list.append(Hop(t*self.time_step*step_skip, structures[t][a+starting_index], 
                                                 a+starting_index, cumul_sum, hopping_length)) 
                                                 # KJ-adding a+starting_index requires that all of Li-site indices are continuous to each other
                                                 # KJ-finsite: the site object when the hop has occured.
                        cumul_sum = np.zeros(3) 

            print('Total number of hops counted: ', len(self.hop_list))
        else:
            self.hop_list = text_input
            print('Imported from text file... Structure and DiffusionAnalyzer objects will not work')


    def export_txt(self, DIR):
        if not len(self.hop_list) == 0:
            with open(DIR, 'w') as f:
                for hop in self.hop_list:
                    f.write('{:0.1f}//{}//{}//{}//{}//{}//{}\n' 
                            .format(hop.time, hop.fin_site.as_dict(), hop.index, 
                                    hop.vector[0], hop.vector[1], hop.vector[2], hop.threshold_length)) 
        else:
            print('No element in hopping list !!')


    def Calculate_HH_Correlation(self, time_scale, length_scale):
        hh_correlation_list = []
#         hh_correlation_number = 0
        for i, ref in enumerate(self.hop_list):
            for j, con in enumerate(self.hop_list):
                if i > j:
                    if ref.index != con.index:
                        if (abs(ref.time-con.time) < time_scale):
#                             if (ref.time < con.time and ref.fin_site.distance(con.ini_site) < length_scale) or \
#                             (ref.time > con.time and ref.ini_site.distance(con.fin_site) < length_scale):
                            if np.average([ref.ini_site.distance(con.ini_site), ref.fin_site.distance(con.fin_site)]) < length_scale:
                                hh_correlation_list.append(HHCorrelation(ref, con, time_scale, length_scale))

                
        return hh_correlation_list

        
    def Calculate_HH_number(self, time_scale, length_scale):
        hh_correlation_number = 0
        for i, ref in enumerate(self.hop_list):
            for j, con in enumerate(self.hop_list):
#                 if i > j:
                    if ref.index != con.index:
                        if (abs(ref.time-con.time) < time_scale):
                            if np.average([ref.ini_site.distance(con.ini_site), ref.fin_site.distance(con.fin_site)]) < length_scale:
                                hh_correlation_number += 1
                                break
        return hh_correlation_number
    
    @classmethod
    def from_txt(cls, DIR, species, temperature):
        try:
            text_hop_list = []
            with open(DIR, 'r') as g:
                for line in g.readlines():
                    t, fs, i, v1, v2, v3, thl = line.split('//')
                    text_hop_list.append(Hop(float(t), PeriodicSite.from_dict(json.loads(fs.replace("'", '"'))), 
                                             int(i), [float(v1), float(v2), float(v3)], float(thl)))
            return cls(None, species, temperature, hopping_length=text_hop_list[0].threshold_length, text_input=text_hop_list)
                    
        except(IndexError):
            print('Text input file is not valid !!')
        

    @classmethod
    def from_paths(cls, paths, species, temperature, step_skip=10, time_step=2, hopping_length=3):
        strs = []
        for x in paths:
            if os.path.exists('{}/vasprun.xml' .format(x)):   
                temp_str = Vasprun('{}/vasprun.xml' .format(x), ionic_step_skip=step_skip, exception_on_bad_xml=False).structures
                strs += temp_str
                del temp_str
                
        return cls(strs, species, temperature, step_skip=step_skip, time_step=time_step, hopping_length=hopping_length)

    
class Rotation:
    def __init__(self, time, duration, angle, index, site, rotation_q, final_q, init_q):
        self.time = time
        self.duration = duration
        self.angle = angle
        self.index = index
        self.site = site
        self.rotation_q = rotation_q
        self.final_q = final_q
        self.init_q = init_q

class RHCorrelation:
    def __init__(self, Rotation, Hop, time_scale, length_scale):
        self.indice = [Rotation.index, Hop.index]
        self.time_scale = time_scale
        self.length_scale = length_scale
        self.timestep = abs(Rotation.time - Hop.time)
        if Rotation.time <= Hop.time:
            self.distance = Rotation.site.distance(Hop.ini_site)
        else:
            self.distance = Rotation.site.distance(Hop.fin_site)
        
    
class RotationAnalyzer:
    def __init__(self, structures, species, temperature, step_skip=10, time_step=2, n_process=None, rot_graph=None, from_init=None, each_time=None, info_dict=None):
        self.species = Element(species)
        self.temperature = temperature
        self.step_skip = step_skip
        self.time_step = time_step
        self.structures = structures
        self.neighbor_sets_matrix = []
        self.rot_graph = []  
        self.centers = []
        self.rotations_from_init = from_init
        self.rotations_each_time = each_time
        self.info_dict = info_dict
        self.split = None
        self.part_index = None

        start_time = 0 
        '''
        start_time : Seems like Byungju made it to continue adding new trajectories. At this moment, disable this.
        '''
        #if rot_graph is not None:
        #    start_time = rot_graph[0, -1, 0] + self.time_step*self.step_skip
        #else:
        #    start_time = 0
        if (not info_dict):
            '''
            Determine Either starting from scratch, or starting from an archived rot graph information
            This case: we have no info_dict, so we start from scratch
            '''
            print("Starting the analysis from scratch")
            # count the number of species
            indices = []
            for i in range(len(structures[0])):
                if structures[0].species[i] == Element(species):
                    indices.append(i)
                    self.centers.append(structures[0][i])
            num_indices = len(indices)
            starting_index = min(indices)
            print('The number of {} atoms: {}' .format(self.species, num_indices))
            # Derive all neighbor_sets objects from the input trajectory       

            with multiprocessing.Pool(processes=n_process) as p:
                self.neighbor_sets_matrix = p.starmap(derive_neighbor_sets,
                                                    [(struct, species, starting_index, num_indices, i, len(structures)) for i, struct in enumerate(structures)])
            self.neighbor_sets_matrix = np.transpose(np.array(self.neighbor_sets_matrix)) #KJ-This is a list of list of 4-coordinate NeighborSet objects
            # KJ-(num_site X num_structure) since it was transposed

#             self.neighbor_sets_matrix = []
#             for i, struct in enumerate(structures):
#                 self.neighbor_sets_matrix.append(derive_neighbor_sets(struct, species, starting_index, num_indices, i, len(structures)))
#             print('Process finished !!')
#             self.neighbor_sets_matrix = np.transpose(np.array(self.neighbor_sets_matrix))

            # Save rotation graph
            '''
            TBD: determine whether we want to reference our quaternions based on t0
            or relaxed position (which we need another Structure object input)
            '''
            for i, tet in enumerate(self.neighbor_sets_matrix): # KJ-For each P-sites
                temp_placeholder = []
                for t, time in enumerate(tet):
                    if time is not None: #Non-zero time
                        #temp_placeholder.append([t*self.time_step*self.step_skip+start_time, calc_reduced_angle(tet[0], time)[0], 
                        #                      calc_reduced_angle(tet[0], time)[1], calc_reduced_angle(tet[0], time)[2]])
                        temp_placeholder.append([t*self.time_step*self.step_skip+start_time, calc_quat(tet[0], time)])
                        #KJ-Appending [tet, [absolute time, angle1, angle2, angle3 in euler angles]]
                    else: #At time=0
                        temp_placeholder.append([t*self.time_step*self.step_skip+start_time, calc_quat(tet[0],tet[0])])
                        
                self.rot_graph.append(temp_placeholder)
            
            if rot_graph is not None:
                self.rot_graph = np.concatenate((rot_graph, np.array(self.rot_graph)), axis=1)
            else:
                self.rot_graph = np.array(self.rot_graph)

            print('Rotation information successfully archived')
            self.info_dict = {}
            #self.info_dict['structures']=self.structures
            self.info_dict['time_step']=self.time_step
            self.info_dict['step_skip']=self.step_skip
            self.info_dict['temperature']=self.temperature
            self.info_dict['species']=self.species
            self.info_dict['centers']=self.centers
            self.info_dict['neighbor_sets_matrix']=self.neighbor_sets_matrix
            self.info_dict['split'] = self.split
            self.info_dict['part_index'] = self.part_index
        elif info_dict['part_index']=='rot_graph_read':
            '''
            Only rotation graph is exported (with info_dict), not from_init or each_time
            '''
            self.info_dict = info_dict
            self.structures = structures
            #self.structures = info_dict['structures']
            self.time_step = info_dict['time_step']
            self.step_skip = info_dict['step_skip']
            self.temperature = info_dict['temperature']
            self.species = info_dict['species']
            self.centers = info_dict['centers']
            self.neighbor_sets_matrix = info_dict['neighbor_sets_matrix']
            self.split = info_dict['split']
            self.part_index = info_dict['part_index']
            self.rot_graph = rot_graph
        else:
            print('Plugging in the pre-analyzed information')
            self.info_dict = info_dict
            self.structures = structures
            #self.structures = info_dict['structures']
            self.time_step = info_dict['time_step']
            self.step_skip = info_dict['step_skip']
            self.temperature = info_dict['temperature']
            self.species = info_dict['species']
            self.centers = info_dict['centers']
            self.neighbor_sets_matrix = info_dict['neighbor_sets_matrix']
            if "split" in info_dict.keys():
                self.split = info_dict['split']
            else:
                self.split = None
            if "part_index" in info_dict.keys():
                self.part_index = info_dict['part_index']
            else:
                self.part_index = None
            self.rot_graph = rot_graph
            #self.rot_graph = text_input[0]
            #self.centers =  text_input[1]
            #print('Imported from text file... neighbor_sets_matrix will not work')
            

    def export_txt(self, DIR):
        if not len(self.rot_graph) == 0:
            with open(DIR, 'w') as f:
                f.write(str(self.rot_graph.tolist()))
                f.write('\n')
                f.write(str([x.as_dict() for x in self.centers]))
        else:
            print('No element in rotation list !!')


    def export_rotation_analysis(self, DIR):

        with open('{}/out_{}.pkl'.format(DIR,"info_dict"), 'wb') as f:
            pickle.dump(self.info_dict, f)
        with open('{}/out_{}.pkl'.format(DIR, "structures"), 'wb') as f:
            pickle.dump(self.structures, f)
        if self.part_index==None:
            if not len(self.rotations_each_time) == 0:
                np.save("{}/out_{}.npy".format(DIR,"each_time"), np.array(self.rotations_each_time), allow_pickle=True)
            if not len(self.rotations_from_init) == 0:
                np.save("{}/out_{}.npy".format(DIR,"from_init"), np.array(self.rotations_from_init), allow_pickle=True)
            #if not len(self.rot_graph) == 0:
            #    np.save("{}/out_{}.npy".format(DIR,"rot_graph"), self.rot_graph, allow_pickle=True) 
        if self.part_index != None:
            if len(self.rotations_each_time) != 0:
                np.save("{}/out_{}_{}.npy".format(DIR,"each_time", np.char.zfill(str(self.part_index),3)), np.array(self.rotations_each_time), allow_pickle=True)
            if len(self.rotations_from_init) != 0:
                np.save("{}/out_{}_{}.npy".format(DIR,"from_init", np.char.zfill(str(self.part_index),3)), np.array(self.rotations_from_init), allow_pickle=True)
            #if len(self.rot_graph) != 0:
            #    np.save("{}/out_{}_{}.npy".format(DIR,"rot_graph", np.char.zfill(str(self.part_index),3)), self.rot_graph, allow_pickle=True) 

    def export_rot_graph_only(self, DIR):
        with open('{}/out_{}.pkl'.format(DIR,"info_dict"), 'wb') as f:
            pickle.dump(self.info_dict, f)
        if not len(self.rot_graph) == 0:
            np.save("{}/out_{}.npy".format(DIR,"rot_graph"), self.rot_graph, allow_pickle=True) 
        
            
    def count_rot_from_graph_from_init(self, max_time_delay=3000, n_process=None):
        rotations=[]
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)
        
        for a, atom in enumerate(self.rot_graph):
            atom_info = []
            with multiprocessing.Pool(processes=n_process) as p:
                atom_info = p.starmap(self.angle_computer_init, [(index, a, num_delayed_frames) for index in range(len(atom)-num_delayed_frames)])
                rotations.append(atom_info)
        self.rotations_from_init = rotations
        return rotations

        '''
        for a, atom in enumerate(self.rot_graph):
            atom_info = []
            with multiprocessing.Pool(processes=n_process) as p:
                time_info = p.starmap(lambdafunc,
                        [])
            for frameindex, frame in enumerate(atom):
                if frameindex >= len(atom)-num_delayed_frames:
                    break
                time_info = []

                for i in range(num_delayed_frames):
                    time_info.append(Quaternion.sym_distance(atom[0][1], atom[frameindex+i][1]))
                atom_info.append(time_info)
            rotations.append(atom_info)
        self.rotations_from_init = rotations
        return rotations
        '''

    def count_rot_from_graph_each_time(self, max_time_delay=3000,n_process=None):
        '''
        outdated, use atomwise class methods
        max_time_delay: dt value maximum (in fs)
        n_process: number of cpu-cores
        '''
        rotations = []
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)

        for a, atom in enumerate(self.rot_graph):
            '''
            Note, rot_graph[P_index][0] = time in fs, rot_graph[1] = Quaternion object (referenced to the start of the simulation
            '''
            atom_info = []
            for frameindex, frame in enumerate(atom):
                if frameindex >= len(atom)-num_delayed_frames:
                    break
                with multiprocessing.Pool(processes=n_process) as p:
                    # atom[0][1] is this P's first snapshot's quaternion object
                    time_info = p.starmap(angle_quat, [(atom[0][1], atom[frameindex+i][1]) for i in range(num_delayed_frames)])
                atom_info.append(time_info)
            rotations.append(atom_info)
                #with multiprocessing.Pool(processes=n_process) as p:
                #    atom_info = p.starmap(self.angle_computer_each_time, [(index, a, num_delayed_frames) for index in range(len(atom)-num_delayed_frames)])
                #    rotations.append(atom_info)
        self.rotations_from_init = rotations
        return rotations

    def count_rot_from_graph_from_init_atomwise(self, max_time_delay=5000,n_process=None, split=1000, part_index=None):
        '''
        max_time_delay: dt value maximum (in fs)
        n_process: number of cpu-cores
        split: for each multiprocessing rounds, how many frames you want to split. i.e. 1000 steps = 10 step_skip*2 time_step = 20 ps
        rotations object: [P-index, frames]
        Memory to handle within the multiprocessing becomes excessive if you make the Pool handle too many trajectories.
        Therefore, we split into "split" number of trajectories each time to analyze.
        '''
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)
        rotation_list = []
        start_index = 0
        split_index = 0
        self.split = split
        self.part_index = part_index
        with multiprocessing.Pool(processes=n_process) as p:
            while (start_index < len(self.rot_graph[0])-num_delayed_frames):
                if (part_index == None) or (part_index == split_index):
                    print('Starmap launching now for from_init analysis')
                    rotation_add = p.starmap(self.lambda_from_init_split, [(i, num_delayed_frames, start_index, split) for i in range(len(self.rot_graph))])
                    rotation_list.append(np.array(rotation_add))
                if split_index == part_index:
                    break
                start_index += split
                split_index += 1
                #p.close()
        #for i in rotation_list:
        #    print(i.shape)
        if len(rotation_list)==0:
            self.rotations_from_init=np.array([])
            return []
        else:
            rotations = np.concatenate([i for i in rotation_list], axis=1)
            #print(rotations.shape)
            self.rotations_from_init = rotations
            return rotations

    def count_rot_from_graph_each_time_atomwise(self, max_time_delay=5000,n_process=None, split=1000, part_index=None):
        '''
        max_time_delay: dt value maximum (in fs)
        n_process: number of cpu-cores
        split: for each multiprocessing rounds, how many frames you want to split. i.e. 1000 steps = 10 step_skip*2 time_step = 20 ps
        rotations object: [P-index, frames]
        Memory to handle within the multiprocessing becomes excessive if you make the Pool handle too many trajectories.
        Therefore, we split into "split" number of trajectories each time to analyze.
        '''
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)
        rotation_list = []
        start_index = 0
        split_index = 0
        self.split = split
        self.part_index = part_index
        with multiprocessing.Pool(processes=n_process) as p:
            while (start_index < len(self.rot_graph[0])-num_delayed_frames):
                if (part_index == None) or (part_index == split_index):
                    print('Starmap launching now for each_time analysis')
                    rotation_add = p.starmap(self.lambda_each_time_split, [(i, num_delayed_frames, start_index, split) for i in range(len(self.rot_graph))])
                    rotation_list.append(np.array(rotation_add))
                if split_index == part_index:
                    break
                start_index += split
                split_index += 1
                #p.close()

        #for i in rotation_list:
        #    print(i.shape)
        if len(rotation_list)==0:
            self.rotations_each_time=np.array([])
            return []
        else:
            rotations = np.concatenate([i for i in rotation_list], axis=1)
            #print(rotations.shape)
            self.rotations_each_time = rotations
            return rotations



    def lambda_from_init_split(self, index, num_delayed_frames, start_index, split):
        '''
        helper function for count_rot_graph_from_init_atomwise
        '''
        atom_info = []
        for frameindex, frame in enumerate(self.rot_graph[index]):
            if frameindex<start_index: continue
                # Setting the start condition for this split analysis
            if frameindex >= len(self.rot_graph[0])-num_delayed_frames:
                break
                # If it is the end of the entire trajectory, finish it here!
            if frameindex >= start_index+split:
                break
                # End of the split, end here
            time_info = []
            for i in range(num_delayed_frames):
                time_info.append(abs_angle_quat(self.rot_graph[index][0][1], self.rot_graph[index][frameindex+i][1]))
            atom_info.append(time_info)
        return atom_info


    def lambda_each_time_split(self, index, num_delayed_frames, start_index, split):
        '''
        helper function for count_rot_graph_from_init_atomwise
        '''
        atom_info = []
        for frameindex, frame in enumerate(self.rot_graph[index]):
            if frameindex<start_index: continue
                # Setting the start condition for this split analysis
            if frameindex >= len(self.rot_graph[0])-num_delayed_frames:
                break
                # If it is the end of the entire trajectory, finish it here!
            if frameindex >= start_index+split:
                break
                # End of the split, end here
            time_info = []
            for i in range(num_delayed_frames):
                time_info.append(abs_angle_quat(self.rot_graph[index][frameindex][1], self.rot_graph[index][frameindex+i][1]))
            atom_info.append(time_info)
        return atom_info

    """
    def count_rot_from_graph_each_time_atomwise(self, max_time_delay=5000,n_process=None):
        '''
        max_time_delay: dt value maximum (in fs)
        n_process: number of cpu-cores
        '''
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)
        with multiprocessing.Pool(processes=n_process) as p:
            rotations = p.starmap(self.lambda_each_time, [(i, num_delayed_frames) for i in range(len(self.rot_graph))])
        self.rotations_each_time = rotations
        return rotations

    def count_rot_from_graph_from_init_atomwise(self, max_time_delay=5000,n_process=None):
        '''
        max_time_delay: dt value maximum (in fs)
        n_process: number of cpu-cores
        '''
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)
        with multiprocessing.Pool(processes=n_process) as p:
            rotations = p.starmap(self.lambda_from_init, [(i, num_delayed_frames) for i in range(len(self.rot_graph))])
        self.rotations_from_init = rotations
        return rotations
        """
    def lambda_from_init(self, index, num_delayed_frames):
        '''
        helper function for count_rot_graph_from_init_atomwise
        '''
        atom_info = []
        for frameindex, frame in enumerate(self.rot_graph[index]):
            if frameindex >= len(self.rot_graph[0])-num_delayed_frames:
                break
            time_info = []
            for i in range(num_delayed_frames):
                time_info.append(abs_angle_quat(self.rot_graph[index][0][1], self.rot_graph[index][frameindex+i][1]))
            atom_info.append(time_info)
        return atom_info

    def lambda_each_time(self, index, num_delayed_frames):
        '''  
        helper function for count_rot_graph_each_time_atomwise
        '''
        atom_info = []
        for frameindex, frame in enumerate(self.rot_graph[index]):
            if frameindex >= len(self.rot_graph[0])-num_delayed_frames:
                break
            time_info = []
            for i in range(num_delayed_frames):
                time_info.append(abs_angle_quat(self.rot_graph[index][frameindex][1], self.rot_graph[index][frameindex+i][1]))
            atom_info.append(time_info)
        return atom_info

    def count_rotations_from_each_time(self, ignored_angles=15, firstmax_detector=15, min_dt=100):
        '''
        parameter ignored_angles=15 : below this value are considered vibrations and smoothed out.
        Set this value to zero if you want no smoothing.
        parameter min_dt=100: the first dt scan, in fs unit (default value: 0.1 fs corresponding to a single vibration of rigid 10E13 attempt freq)
        Algorithm:
        1. generate a t0 vs dt=100 plot (y-axis=angles in degrees)
        2. Make a list of the t0, which is the onset of each peaks with angle larger than ignored_angles
        3. For each t0, keep increasing dt until the value sees a drop and find the maximum angle at which angle starts to decrease. 
            For this task, use scipy.signal.find_peaks with respect to dt, then obtain the information about the first peak.
        4. (Optional): peak['peak_heights'], but this gives you the center of the peak position.
        Note, we can also use "scipy.find_peaks_cwt"
        Note, we can also use "scipy.signal.peak_widths"
        Note, now I think that using "ruptures" package to detect change is better.
        If this does not work, set a very rough and simple cutoff.
        Since ruptures can work with an unspecified number of changes, it is not 100% reliable and every time
        the parameters may have to be checked. Therefore, let's use a very crude method.
        parameter ; threshold. dt keep increasing, record the maximum value (cumul_max_dt) at each dt
        until you meet a dt where angle < cumul_max_dt-threshold. Then, return cumul_max_dt and the time when it occured.
        TBD: implement this!
        '''
        indices = []
        for i in range(len(self.structures[0])):
            if self.structures[0].species[i] == Element(self.species):
                indices.append(i)
        num_indices = len(indices)
        starting_index = min(indices)

        num_min_dt_frames = min_dt/self.step_skip/self.time_step
        assert num_min_dt_frames.is_integer()
        num_min_dt_frames = int(num_min_dt_frames)
        full_init_scan = []
        full_init_peak = []
        for P_index, P in enumerate(self.rotations_each_time):
            init_scan = []
            init_peak = []
            for t0_index, t0_frame in enumerate(P):
                val = t0_frame[num_min_dt_frames]
                val = val*180/np.pi
                # CONVERTING FROM RADIAN TO DEGREES
                if val<ignored_angles:
                    val = 0
                init_scan.append(val)
            onset = False
            for t0_index, t0_angle in enumerate(init_scan):
                if (t0_angle>0) and (onset == False):
                    # Start of a rotation!
                    # FYI, continuation of a rotation would be t0_angle >0 and onset = True
                    onset = True
                    init_peak.append(t0_index)
                    continue
                if (t0_angle==0) and (onset == True):
                    onset = False
                    # Finish of a rotation
            full_init_scan.append(init_scan)
            full_init_peak.append(init_peak)
        # Document the position of initial peaks.
        # Now, for each peak position, we change dt and find the maximum value here.
        rot_out = []
        for P_index, P in enumerate(full_init_peak):
            P_rots = []
            for init_peak in P:
                interested_t0 = self.rotations_each_time[P_index][init_peak] # In 1D list of radians with increasing dt at fixed t0
                firstmaxangle, firstmaxdt = self.maxangle(interested_t0, firstmax_detector) # In degrees!
                # Now that we know when the rotation occured, let's record this rotation event
                # as a Rotation object.
                event_frame = init_peak+num_min_dt_frames
                duration_frame = firstmaxdt-num_min_dt_frames
                event_time = event_frame*self.step_skip*self.time_step # in fs
                # duration is recorded as the time to reach the maximum angle
                duration = duration_frame*self.step_skip*self.time_step # in fs
                final_q = self.rot_graph[P_index][event_frame+duration_frame][1]
                init_q = self.rot_graph[P_index][event_frame][1]
                rot_quat = final_q * init_q.inverse
                # Q(final) = P(diff) * I(init)
                # P(diff) = Q(final)*I(init).inverse
                # we have [1] because rot_graph[P_index][time_index] = np.array([time in fs, Quaternion object])
                r = Rotation(time=event_time, duration=duration, angle=firstmaxangle, index=P_index+starting_index,
                             site=self.structures[init_peak+num_min_dt_frames][P_index+starting_index],
                             rotation_q=rot_quat, final_q=final_q, init_q=init_q)
                P_rots.append(r)
            rot_out.append(P_rots)
        return rot_out


    def maxangle(self, ar, tolerance=15):
        '''
        ar: 1D array of angles at t0, t0+dt, varying dt
        tolerance: decreasing angle tolerance in degrees
        this method finds the first maximum peak before any decrease is detected within the tolerance angle

        Returns (max_so_far, max_index)
        max_so_far: the maximum angle of rotation in degrees
        max_index: the time it took to reach the maximum angle of rotation
        '''
        ar_in_deg = ar.copy()
        ar_in_deg = ar_in_deg*180/np.pi
        max_so_far = ar_in_deg[0]
        max_index = 0
        for index, i in enumerate(ar_in_deg):
            if i>max_so_far:
                max_so_far = i
                max_index = index
            if i<max_so_far-tolerance:
                return (max_so_far, max_index)
        return (max_so_far, max_index)

    '''
    IDEA; once we compute the hops and rotations and find that certain hops occur within temporal and spatial proximity to a rotation event, let's try to debunk the notion
    that PS4 pushes the lithium. 
    Based on the center P atom, assuming that the lithium exactly follows the PS4 rotation, we can find the hypothetical coordinate of the Li_rotated.
    Li_rotated - Li_init vs. Li_hopped_Li_init -> For these two vectors, we can compute the angle between these two vectors. If angle is closer to zero, then there is a correlation between the PS4
    rotation direction and Li-hop.        
    '''


    '''
        for a, atom in enumerate(self.rot_graph):
            atom_info = []
            for frameindex, frame in enumerate(atom):
                if frameindex >= len(atom)-num_delayed_frames:
                    break
                with multiprocessing.Pool(processes=n_process) as p:
                    time_info = p.starmap(angle_quat, [(atom[0][1], atom[frameindex+i][1]) for i in range(num_delayed_frames)])
                atom_info.append(time_info)
            rotations.append(atom_info)
                #with multiprocessing.Pool(processes=n_process) as p:
                #    atom_info = p.starmap(self.angle_computer_each_time, [(index, a, num_delayed_frames) for index in range(len(atom)-num_delayed_frames)])
                #    rotations.append(atom_info)
        self.rotations_from_init = rotations
        return rotations
       
        rotations = []
        num_delayed_frames = int(max_time_delay/self.step_skip/self.time_step)

        for a, atom in enumerate(self.rot_graph):
            atom_info = []
            for frameindex, frame in enumerate(atom):
                if frameindex >= len(atom)-num_delayed_frames:
                    break
                time_info = []
                for i in range(num_delayed_frames):
                    time_info.append(Quaternion.sym_distance(frame[1], atom[frameindex+i][1]))
                atom_info.append(time_info)
            rotations.append(atom_info)
        self.rotations_each_time = rotations
        return rotations
    '''

    def plot_rotations(self, mode, P_index='all',vmax=180, filename=None):
        '''
        Plot t0-dt diagram
        Functionalities to add:
        (1) split into 50 ps diagrams (250 ps simulation would end up with 5 * 50 ps diagrams)
        (2) if we know the time when rotation was detected, also plot those information here
        '''
        if mode=='from_init':
            rotations = self.rotations_from_init
        elif mode=='each_time':
            rotations = self.rotations_each_time
        else:
            print('Wrong mode, please double check')
        for rot_index, rot in enumerate(rotations):
            if not ((P_index =='all') or (P_index==rot_index)):
                continue
            print("Plotting {}th-index PS4".format(rot_index))
            x = []
            y = []
            z = []
            for index_i, i in enumerate(rot):
                for index_j, j in enumerate(i):
                    x.append(index_i*self.step_skip*self.time_step/1000)
                    y.append(index_j*self.step_skip*self.time_step/1000)
                    z.append(j*180/np.pi)
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            xy = np.column_stack([x.flat, y.flat])
            print(xy)
            xmin, xmax, ymin, ymax = 0,max(x),0,max(y)
            grid_x, grid_y = np.mgrid[xmin:xmax:1000j, ymin:ymax:1000j]
            method = 'linear'
            fig, ax = plt.subplots(figsize=(30, 30*max(y)/max(x)))
            grid_z = scipy.interpolate.griddata(xy,z,(grid_x, grid_y), method=method, fill_value=0.2)
            # [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
            import numpy.ma as ma
            im = ax.pcolormesh(grid_x, grid_y, ma.masked_invalid(grid_z), cmap='gist_earth', vmin=0, vmax=vmax)
            #contours=[10,20,30,40,50,60]
            #contour_z = contours
            #ax.contour(grid_x, grid_y, grid_z, levels=contour_z, linewidths=4, colors='black', linestyles='dashed')
            #im = ax.pcolormesh(grid_x, grid_y, ma.masked_invalid(grid_z), cmap='jet', vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z))
            fig.colorbar(im, ax=ax)
            ax.tick_params(axis='both',  reset=True, which='both', direction='out', length=10, width=3, color='black', top=False, right=False, zorder=100000)
            ax.set_xlabel(r'$t_0$ (ps)')
            ax.set_ylabel(r'$dt$ (ps)')
            if filename:
                fig.savefig(filename+'.pdf', dpi=300)
            return
        """
        angle = np.radians(self.rot_graph[:, 1:, 1:]) # Indices: P-site, time, angle_information
        # Last index is [1:] since the 0-th component is time
        disp = np.radians(self.rot_graph[:, 1:, 1:] - self.rot_graph[:, :-1, 1:]) 
        # disp is subtracting by a single timestep, so angular velocity at each moment
        disp = (disp/np.pi - np.round(disp/np.pi))*np.pi # radian / fs
        self.disp = disp

        self.simulation_time = len(structures)*self.step_skip*self.time_step
        '''
        conv_disp = np.stack([[(disp[i, :, 2]*np.sin(angle[i, :, 0])*np.sin(disp[i, :, 1]) + disp[i, :, 0]*np.cos(disp[i, :, 1])) for i in range(angle.shape[0])],
             [(disp[i, :, 2]*np.sin(angle[i, :, 0])*np.sin(angle[i, :, 1]) - disp[i, :, 0]*np.cos(angle[i, :, 1])) for i in range(angle.shape[0])],
             [(disp[i, :, 2]*np.cos(angle[i, :, 0]) + disp[i, :, 1]) for i in range(angle.shape[0])]], axis=-1)
        x: term1 and 0 flipped
        y: cos was wrongly written as sin, sign is opposite
        z: term 1 and 0 filpped
        '''
        conv_disp = np.stack([[(disp[i,:,2]*np.sin(angle[i,:,1])*np.sin(disp[i,:,0]) + disp[i,:,1]*np.cos(disp[i,:,0])) for i in range(angle.shape[0])],
                            [(-disp[i,:,2]*np.sin(angle[i,:,1])*np.cos(angle[i,:,0])+disp[i,:,1]*np.sin(disp[i,:,0])) for i in range(angle.shape[0])],
                            [(disp[i,:,2]*np.cos(angle[i,:,1])+disp[i,:,0]) for i in range(angle.shape[0])]], axis=-1)
        # angular velocity equation: https://www.lehman.edu/faculty/dgaranin/Mechanics/Mechanis_of_rigid_bodies.pdf
        # disp is velocity values. conv_disp is converted angular velocity projected to ex, ey, ez axis.
        indices = []
        for i in range(len(structures[0])):
            if structures[0].species[i] == Element(self.species):
                indices.append(i)
        starting_index = min(indices)        

        self.rotation_list = []

        for a, atom in enumerate(conv_disp):
            cumul_sum = np.zeros(3)
            t_prev = 0
            for t, time in enumerate(atom):
                cumul_sum = cumul_sum + time
                if np.linalg.norm(cumul_sum) > rotation_threshold:
                    self.rotation_list.append(Rotation(t*self.time_step*self.step_skip, (t-t_prev)*self.time_step*self.step_skip, structures[t][a+starting_index], 
                                             angle[a, t], a+starting_index, cumul_sum, rotation_threshold))
                    cumul_sum = np.zeros(3)
                    t_prev = t
        return len(self.rotation_list)
        """

          
    def Calculate_RH_Correlation(self, hop_list, time_scale, length_scale):
        rh_correlation_list = []
        for i, hop in enumerate(hop_list):
            for j, rot in enumerate(self.rotation_list):
                if abs(hop.time-rot.time) <= time_scale:
                    if (hop.fin_site.distance(rot.site) <= length_scale) or (hop.ini_site.distance(rot.site) <= length_scale):
                        rh_correlation_list.append(RHCorrelation(rot, hop, time_scale, length_scale))
          
        return rh_correlation_list

#     def Calculate_RH_Uncorrelation(self, hop_list, uncorr_time=50000, uncorr_length=6):
# #         hh_uncorrelation_number = len(self.hop_list)
#         rh_uncorrelation_number = 0
        
#         for i, hop in enumerate(hop_list):
#             for j, rot in enumerate(self.rotation_list):
#                 if (i > j and (abs(hop.time-rot.time) > uncorr_time) and \
#                 (hop.fin_site.distance(rot.site) > uncorr_length) and (hop.ini_site.distance(rot.site) > uncorr_length)):
# #                     hh_uncorrelation_number -= 1
# #                     break
#                     hh_uncorrelation_number += 1
    
    def Calculate_RH_number(self, hop_list, time_scale, length_scale):
        rh_correlation_number = 0
        for hop in hop_list:
            for rot in self.rotation_list:
                if abs(hop.time-rot.time) < time_scale:
                    if np.average([hop.fin_site.distance(rot.site), hop.ini_site.distance(rot.site)]) <= length_scale:
                        rh_correlation_number += 1
                        break
                
        return rh_correlation_number

        
#     def Calculate_RH_Uncorrelation(self, hop_list, uncorr_time=5000, uncorr_length=6):
#         rh_uncorrelation_number = 0
 
#         for hop in hop_list:
#             counter = 0
#             for rot in self.rotation_list:
#                 if abs(hop.time-rot.time) < uncorr_time:
#                     if ((hop.time <= rot.time) and (hop.fin_site.distance(rot.site) <= uncorr_length)) or \
#                     ((hop.time > rot.time) and (hop.ini_site.distance(rot.site) <= uncorr_length)):
#                         counter += 1
#             if counter == 0:
#                 rh_uncorrelation_number += 1
                
#         return rh_uncorrelation_number    
    def angle_computer_init(self, t0_index, atom_index, num_delayed_frames):
        time_info = []
        for i in range(num_delayed_frames):
            time_info.append(Quaternion.sym_distance(self.rot_graph[atom_index][0][1], self.rot_graph[atom_index][t0_index+i][1]))
        return time_info
    def angle_computer_each_time(self, t0_index, atom_index, num_delayed_frames):
        time_info = []
        for i in range(num_delayed_frames):
            time_info.append(Quaternion.sym_distance(self.rot_graph[atom_index][t0_index][1], self.rot_graph[atom_index][t0_index+i][1]))
        return time_info
    
    @classmethod
    def from_txt(cls, DIR, species, temperature):
        try:
            text_rot_list = []
            with open(DIR, 'r') as g:
                lines = g.readlines()
                text_rot_list = json.loads(lines[0].replace('nan', 'NaN'))
                text_center_list = json.loads(lines[1].replace("'", '"'))

            return cls(None, species, temperature, text_input=[np.array(text_rot_list), [PeriodicSite.from_dict(x) for x in text_center_list]])
                    
        except(IndexError):
            print('Text input file is not valid !!')

    @classmethod
    def from_npys(cls, DIR):
        print('Reading single npy file')
        rot_graph = np.load("{}/out_{}.npy".format(DIR, "rot_graph"), allow_pickle=True)
        from_init = np.load("{}/out_{}.npy".format(DIR, "from_init"), allow_pickle=True)
        each_time = np.load("{}/out_{}.npy".format(DIR, "each_time"), allow_pickle=True)
        with open('{}/out_{}.pkl'.format(DIR, "info_dict"), 'rb') as f:
            info_dict = pickle.load(f)
        return cls(None, species=info_dict['species'], temperature=info_dict['temperature'],info_dict=info_dict, from_init=from_init,each_time=each_time, rot_graph=rot_graph)

    @classmethod
    def from_many_npys(cls, DIR, read_structures=False):
        print('Reading multiple npy files and concatenating them')
        '''
        Read multiple npy files to make a combined class!
        '''
        def read_output(type_output, DIR):
            '''
            Helper function to read multiple split output files (npy files)
            '''
            path = [x for x in os.listdir(DIR) if x.startswith('out_{}'.format(type_output))]
            path.sort()
            full_paths = ['{}/{}' .format(DIR, x) for x in path]
            result_lists = [np.load(i, allow_pickle=True) for i in full_paths]
            result_combined = np.concatenate(result_lists,axis=1)
            return result_combined
        print("Reading rot_graph")
        rot_graph = read_output('rot_graph',DIR)
        print("Reading from_init")
        from_init = read_output('from_init', DIR)
        print("Reading each_time")
        each_time = read_output('each_time', DIR)
        print("Reading info_dict")
        with open('{}/out_{}.pkl'.format(DIR, "info_dict"), 'rb') as f:
            info_dict = pickle.load(f)
        info_dict['part_index']='combined'
        if read_structures:
            print("Reading structures")
            with open('{}/out_{}.pkl'.format(DIR, "info_dict"), 'rb') as f:
                structures = pickle.load(f)
        else:
            structures = None
        print("RotationAnalyzer ready")
        return cls(structures, species=info_dict['species'], temperature=info_dict['temperature'],info_dict=info_dict, from_init=from_init,each_time=each_time, rot_graph=rot_graph)

    @classmethod
    def from_rot_graph(cls, DIR):
        print('Reading pre-analyzed rot_graph')
        rot_graph = np.load("{}/out_{}.npy".format(DIR, "rot_graph"), allow_pickle=True)
        with open('{}/out_{}.pkl'.format(DIR, "info_dict"), 'rb') as f:
            info_dict = pickle.load(f)
        info_dict['part_index']='rot_graph_read'
        return cls(None, species=info_dict['species'], temperature=info_dict['temperature'],info_dict=info_dict, from_init=None,each_time=None, rot_graph=rot_graph)

    
    @classmethod
    def from_paths(cls, paths, species, temperature, step_skip=10, time_step=2, n_process=None, rot_graph=None):
        print('Reading structures from list of paths')
        strs = []
        for x in paths:
            if os.path.exists('{}/vasprun.xml' .format(x)):   
                temp_str = Vasprun('{}/vasprun.xml' .format(x), ionic_step_skip=step_skip, exception_on_bad_xml=False, parse_potcar_file=False).structures
                strs += temp_str
                del temp_str
        
        return cls(strs, species, temperature, step_skip=step_skip, time_step=time_step, n_process=n_process, rot_graph=rot_graph)    



# %%
