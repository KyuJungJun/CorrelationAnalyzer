from pymatgen.io.vasp.outputs import Oszicar
import os

if __name__=="__main__":
    DIR = '.'
    paths = [x for x in os.listdir(DIR) if x.startswith('run')]
    paths.sort()
    full_paths = ['{}/{}' .format(DIR, x) for x in paths]
    print('Reading all Oszicars in this folder')
    num_ionic_steps = 0
    for x in paths:
        if os.path.exists('{}/OSZICAR'.format(x)):
            osz = Oszicar('{}/OSZICAR'.format(x))
            num_ionic_steps += len(osz.ionic_steps)
    print('Total ionic steps: {} steps'.format(num_ionic_steps))
