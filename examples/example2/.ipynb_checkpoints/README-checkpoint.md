Example 2:
CO13 Smart Excess Template Binned Likelihood fit

Description:
This example does a single binned likelihood fit of just the hydrogen gas component, comparing how likely
gamma rays from hydrogen gas (H2) distributions modeled from machine-learning dervied CO13 measurements predict 
our observed data compared to gamma rays from hydrogen gas (H2) distributions modeled from CO12 measurements. 

We use a convolutional neural network defined in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.063018
to predict the CO13 abundance (observed by MOPRA) from the CO12 distribution (observed by DAME) for a small
1deg x 50deg region of the sky along the galactic plane.

Ultimately, we want to use this to predict the CO13 abundance near the galactic center, where we only have 
CO12 abundance data.

There are many sources of gamma rays in the sky, one component comes from pions interacting with hydrogen 
ga (H2). To construct an accurate model of the gamma rays from the H2 component, we need an accurate model 
of the column density of H2 along our line-of-sight (LOS). Since hydrogen gas is a very stable molecule and 
does not emit its own light by itself, we must infer the H2 density by using a proxy gas (carbon monoxide 
CO12), that exists in a proportional abundance (e.g. rho(H2) = rho(CO12) * const. factor. )

However, this gas has some density and when they emit light, some of it can get reabosorbed in optically 
thick clouds along the LOS, and we end up undercounting the amount of light. This specifically biases the 
densest regions, such as in places where we really want to know the distribution of dark matter! So we 
really want to use a less dense gas such as another isotopologue of CO12 (CO13 or CO16/18, but sometimes 
CO16/CO18 are so much more sparse that is is hard to accurately measure their abundances, so we will stick 
with CO13 in this example).

We have maps of CO13 from the MOPRA space telescope and maps of CO12 from the DAME space telescope.
These maps were then used to predict the H2 maps, and then we predicted the gamma-rays we would have observed 
from that component using a software called GALPROP. This gives us maps of the CO13, CO12 derived components 
of the gamma-ray sky, broken down into 17 galactocentric annulus bins (0-16) and 21 different energy bins.

Prior to this code, we have bined these maps into 4 annuli (0-5, 6-9,10-12,13-16) and integrated over all 
energies. These maps are stored under the /data directory.

*The goal of this example is to understand how much more significant are the excess gamma-rays produced by using 
CO13 as a predictor for the H2 component compared to using CO12 as a predictor for the H2 component.*

Here we construct a null-hypoithesis test.

We try to fit the gamma-ray maps predicted from CO12 to the CO13 predicted maps (by varying the normalization of 
the maps at each annulus); this is the "null" model. Then we try to fit the gamma-ray maps of CO12 + a single 
integrated image of the excess gamma rays predicted by using CO13 (the "modified excess template") to the same CO13 
predicted maps; this is the "alternative" model. And we find the 
relative significance using of the log-likelihood. The code here will do 1 iteration of this, but generally we would 
do 1000 realizations (by allowing the CO13 predicted maps that we fit to to vary with noise at each iteration.)

The output of each iteration is the "test-statistic" or "TS". 
TS = -2*(null - alternative)

The sqrt(TS) tells us the sigma-significance of our excess gamma-rays above the gamma-rays prediceted by using 
the canonical model from CO12.

Instructions:

To run this script, you can simply do:

(1) python run.py <sim #>

(2) Run from multitest.sub to run it on the compute nodes of this cluster (preferable).

(3) If you want to run many sims, run from run_jobs.py which will call multitest.sub to run 1000 sims in batches 
    of 50 at a time so that the HPC people don't get mad at you for hogging all the cluster CPUS :/
    (speaking from experience, oops)

Tips:
Fermipyruns many tasks in each simulation and generates many intermediate files, some of which can take up quite
a bit of space on the system. If you were also to run each sim from start to finish each time, each sim would take 
~6hrs to run too. These intermediate files define the config of the simulation, the exposure map/response function
of the telescope, the current data, etc. The good news is that some of these files can be re-used between different 
simulations.

The fermipy package will check if these files exist, and will only run the intermediate tasks and re-generate these 
files if they don't already exist in the working directory of the simulation. If you were to run the simulations using 
the saved intermediate files, each simulation will only take ~10min.

I've included some lines of code in here to help with the workflow for that!

My process is this:

(1) Create a work directory (e.g. workdir_example) that has a config.py file and 
    a run.py file. I'll typically follow the naming convention sims_<workdir name> 
    (e.g. sims_example). The output directory should contain a directory called 
    /out. This code will save your sim data and output maps to the corresponding 
    output folder /output/sim_#, and the log files will be under /out/out-#.txt, 
    and the TS data will be under /tsout/TS_sim_#.txt. Set these parameters in 
    this file at the beginning of main().
    
(2) Run one simluation using multitest.sub (this will be sim 0), running a full
    simulation from start to finish (6hrs), saving all the intermediate files
    along the way. Ensure that lines 167-170 are *uncommmented* and line 335 is *commented*.
    
    $ sbatch --array=0-0 multitest.sub

(3) Once that finishes go to the /output/sim_0 directory and copy all the 
    files over to the work directory 
    
    $ cp *.fits ../../workdir_example/
    
    Then copy the ltcube*.fits file and save it to the fermi-ml/data directory.
    
(4) If that worked, yay! Go ahead and ensure that lines 167-170 are *commmented* and 
    line 335 is *uncommented*. The run 1000 sims by doing:
    
    $ python run_jobs.py
    
    (It is good practice to run 1 sim first to make sure e
    everything works properly before running a thousand or so
    anyways).