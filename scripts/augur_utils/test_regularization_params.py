import os
import sys
import subprocess
import colored_traceback.always
import numpy as np

lam_pot = np.linspace(0., 1.5, 5) # Pick 5 values per parameter, spanning ~150% of the default value on either side
lam_avi = np.linspace(0., 7.5, 5)
lam_drop = np.linspace(0., 2.5, 5)

counter = 1

for p in lam_pot:
    for a in lam_avi:
        for d in lam_drop:
            print (p, a, d, 'Iteration %d of %d'%(counter, len(lam_pot)*len(lam_avi)*len(lam_drop)))
            subprocess.check_call("python dengue/dengue_titers.py -s 2 -v 3 --lam_pot %f --lam_avi %f --lam_drop %f --training_fraction 0.9 --load"%(p, a, d), shell=True)
            counter += 1
