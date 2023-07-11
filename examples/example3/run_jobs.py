import os
import time

batches = [i+1 for i in range(0,1000,50)]
for batch in batches:
    os.system("sbatch --array={}-{} multi-test.sub".format(str(batch), str(batch+50)))
    time.sleep(600)
    os.system("scancel -u abrought")
    os.system("find /pub/abrought/fermi-ml/examples/example3/output/ -type f -name \'ft*.fits\' -delete")
    os.system("ls -l | wc -l")
