import pickle
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
import os
slist = []
for root, dirs, files in os.walk("."):
	filelist = files
	filelist.sort(key=lambda x:int(x.split("_")[1].split(".")[0]))
	for file in filelist:
		print(file)
		with open("{}/{}".format(root, file), "rb") as f:
			vr = pickle.load(f)
			slist.extend(vr.structures)
	break
with open("trajectory.pickle", "wb") as f:
	pickle.dump(slist, f)
print("Generating structure pickle file successful!")