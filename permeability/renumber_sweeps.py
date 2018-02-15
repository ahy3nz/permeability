import numpy as np
import subprocess
import os

all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and 'sweep' in thing[0:5]]
for i, old_folder in enumerate(all_sweeps):
    new_folder = "sweep{}".format(i)
    p = subprocess.Popen("mv {0} {1}".format(old_folder,new_folder), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    #print("moving {} to {}".format(old_folder, new_folder))

