import glob
import re
import numpy as np
import os

simulation="sims_w_gp_final_excess"

files = glob.glob("./" + simulation + "/tsout/*.txt")
TSs = []
sims = []

for f in files:
    nums = [int(s) for s in re.split(r"[_.\s]\s*",f) if s.isdigit()]
    sims.append(nums[0])
    with open(f, 'r') as fp:
        x = [float(line) for line in fp]
        TSs.append(x[0])

TSs = np.asarray(TSs)
sims = np.asarray(sims)

for sim in sims[np.where(TSs < 16.)]:
    print("Removing sim #" + str(sim))
    os.system("rm ./" + simulation + "/tsout/*" + str(sim) + "*.txt")
    os.system("rm -r ./" + simulation + "/output/sim_" + str(sim))
    os.system("rm ./" + simulation + "/*" + str(sim) + "*.dat")

print("Done.")

