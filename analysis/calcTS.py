import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

targetdir = 'sims_MapCubeET_smart_testing_region'
files = sorted(glob.glob("/pub/abrought/h2sim_output/%s/tsout/*.txt" %targetdir))

ts2 = []
for f_name in files:
    f = open(f_name, "r")
    ts2.append(float(f.readline()))

ts = [np.sqrt(val) if np.isfinite(val) else np.nan for val in ts2]

print("mean: " + str(np.nanmean(ts)))
print("median: " + str(np.nanmedian(ts)))
print("N = " + str(len(ts)))

