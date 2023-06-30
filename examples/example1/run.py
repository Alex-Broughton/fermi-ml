"""
Example 1:
CO13  Modified Excess Template Binned Likelihood fit

Description:
This example does a single binned likelihood fit of just the hydrogen gas component, comparing how likely
gamma rays from hydrogen gas (H2) distributions modeled from CO13 measurements predict our observed data 
compared to gamma rays from hydrogen gas (H2) distributions modeled from CO12 measurements. 

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

(3) Once that finishes go to the ~/sims_example/output/sim_0 directory and copy all the 
    files over to the work directory 
    
    $ cp *.fits ../../workdir_example/
    
    Then copy the ltcube*.fits file and save it to the fermi-ml/data directory.
    
(4) If that worked, yay! Go ahead and ensure that lines 167-170 are *commmented* and 
    line 335 is *uncommented*. The run 1000 sims by doing:
    
    $ python run_jobs.py
    
    (It is good practice to run 1 sim first to make sure e
    everything works properly before running a thousand or so
    anyways).

"""


# Imports
from time import sleep
from random import randint
import resource
import random,shutil,yaml
import os,sys
from math import *
import numpy as np
import matplotlib as matplotlib
matplotlib.use('agg')
from fermipy.gtanalysis import GTAnalysis
import pyLikelihood
from BinnedAnalysis import *
import pandas as pd
import glob

def CalcFlux(Fit,name):
    # Helper function
    flux_list = []
    energy_mid_list = []

    energy_list = Fit.energies.tolist()
    for E in range(0,len(energy_list)-1):

        E_low = energy_list[E]
        E_high = energy_list[E+1]
        energy_mid = (E_low + E_high) / 2.0
        flux =  Fit.flux(name,E_low,E_high,energyFlux=True) #MeV/cm^2/s
        flux_list.append(flux)
        energy_mid_list.append(energy_mid)

    #print flux_list
    #print energy_list

    return energy_mid_list,flux_list

def main(cmd_line):
    # Initializations
    
    sim = cmd_line[1]
    indir = os.chdir()
    outdir = "%s/output/sim_%s" % (%indir, %sim)

    # Check that the output directory exists, if not create one.
    # The out put directory should have the format "sims_<name of workdir>"
    # and contain a directory for the output sims, the output logs (out) 
    # and the test statistic outputs (tsout)
    
    if(!os.path.isdir("%s/output" %indir)):
        os.mkdir("mkdir %s/output" %indir)
    if(!os.path.isdir("%s/tsout" %indir)):
        os.mkdir("mkdir %s/tsout" %indir)
    if(!os.path.isdir("%s/spectraout" %indir)):
        os.mkdir("mkdir %s/spectraout" %indir)
    if(os.path.isdir(outdir)==True):
        shutil.rmtree(outdir)
    os.system('mkdir %s' %outdir)
    os.chdir(outdir)

    shutil.copy2('%s/config.yaml' %indir, 'config.yaml') # Always copy config.yaml
    # The sims generate some intermediate files, which define the config of the
    # simulation and the exposure map of the telescope, the data, etc.
    # Some of these files can be re-used between different simulations.
    # The fermipy package will re-generate these files if they do not exist,
    # so we can copy them to the output sim directory so that it won't regenerate
    # them. The first time we run this code, we should generate all of them and 
    # all of these lines should be uncommented (~6hr run time for 1 sim). Subsequently,
    # we can re-use these data products and it will only take (~10 min per sim).
    #shutil.copy2('%s/srcmap_00.fits' %indir, 'srcmap_00.fits')
    #shutil.copy2('%s/bexpmap_00.fits' %indir, 'bexpmap_00.fits')
    #shutil.copy2('%s/ccube_00.fits' %indir, 'ccube_00.fits')
    #shutil.copy2('%s/ft1_00.fits' %indir, 'ft1_00.fits')

    
    
    # SETUP THE ANALYSIS
    
    # Load in the config file, which defines the parameters of our binned likelihood fit
    # This will define the energy range, region of the sky, etc. to fit over.
    # Also in the config file are the maps of all the components we need. They are listed 
    # in a specific order, and numbered as "galdiff##" ("galactic diffuse component ##")
    # Each file corresponds to a map from a specfic source (CO12/CO13) and annulus (0-5, 
    # 6-9,10-12,13-16). We can reference each file by this name. We also have maps ("excess 
    # templates") that are the total integrated excess of gamma-ray (integrated over all 
    # annulus bins) that is predicted by using the better CO13 method over the poor CO12 
    # method.
    
    gta = GTAnalysis('config.yaml',logging={'verbosity' : 3})
    gta.setup()
    #gta.load_roi("after_setup")

    # All data maps are loaded in the config, and can be turned on/off 
    # by fixing there normalization.
    # We will start by turning off all gamma-ray components of the sky 
    # except for the component from atomic (H2) emission, as predicted 
    # from the CO13 distribution rho(H2) = rho(CO13) * const. 
    # The default normalization is just 1 here.
    # This is essentially the "truth" gamma-ray map.
    gta.set_norm("galdiff04",0.0) #CO12_r=0-5
    gta.set_norm("galdiff05",0.0) #CO12_r=6-9
    gta.set_norm("galdiff06",0.0) #CO12_r=10-12
    gta.set_norm("galdiff07",0.0) #CO12_r=13-16
    gta.set_norm("galdiff08",0.0) 
    gta.set_norm("galdiff09",0.0) 
    
    # With just the "truth" maps loaded and non-zero, we will create one 
    # "realization" of Fermi-LAT observations by randomizing ("simulating") 
    # this with a Poisson distribution (assume the counts in each pixel follow 
    # a normal distribution with sigma=sqrt(counts), and randomly sample from 
    # the distribution for each pixel). Each realization will perform a 
    # lieklihood fit on the same maps, but with different noise to simulate
    # a different observation.
    gta.simulate_roi(randomize=True)
    
    #This line will save the "simulated" truth map at this stage
    #gta.write_model_map('modifiedmap_sim_%s' %sim) 
    
    # Delete "truth" sources that were simulated (they are no longer needed).
    #gta.delete_source("MapSource",delete_source_map=False)
    gta.delete_source("galdiff00",delete_source_map=False)
    gta.delete_source("galdiff01",delete_source_map=False)
    gta.delete_source("galdiff02",delete_source_map=False)
    gta.delete_source("galdiff03",delete_source_map=False)

    
    
    # NULL TEST
    
    # Now, we will turn off all gamma-ray components of the sky 
    # except for the component from atomic (H2) emission, as predicted 
    # from the CO12 distribution (the bad way of predicting H2 since CO12 
    # is denser than CO13, and we get more reabsorption along the 
    # line-of-sight (LOS), which results in under predicting H2 in the densest
    # regions of the sky)
    # i.e. rho(H2) = rho(CO12) * const. 
    
    # We can keep the default normalization of 1 to start, or we can start with 
    # some gaussian-random initial normalization.
    
    #set random	 normalizations of sources for performing fit:
    #n4 =  np.random.normal(1.0,0.2)
    #n5 =  np.random.normal(1.0,0.2)
    #n6 =  np.random.normal(1.0,0.2)
    #nms = np.random.normal(1e-4,0.5e-4)
    gta.set_norm("galdiff04",0.8)
    gta.set_norm("galdiff05",0.8)
    gta.set_norm("galdiff06",1.2)
    gta.set_norm("galdiff07",1.0)

    # Now we will fit these gamma ray maps from the H2 component as derived by CO12 
    # to the "truth" dataset by varying the normalization of each map in each annulus.
    
    # So, free all other parameters and fit:
    gta.free_sources(free=True)
    gta.free_source("galdiff07",free=False) # Very little flux in this coming from this annulus, and we ignore it by just fixing this to be constant
    gta.free_source("galdiff08",free=False)
    gta.free_source("galdiff09",free=False)

    Fit = gta.fit()
    null = Fit["loglike"] # this gives us our likelihood for our null test!
    gta.set_norm("galdiff04",0.0) # CO13 r=0-5
    gta.set_norm("galdiff05",0.0) # CO13 r=6-9
    gta.set_norm("galdiff06",0.0) # CO13 r=10-12
    gta.set_norm("galdiff07",0.0) # CO13 r=13-16
    gta.set_norm("galdiff08",0.0) # CO13 excess template map
    gta.set_norm("galdiff09",0.0) # Smart excess template map
    # Save the best-fit map to a file.
    gta.write_model_map('null_model_sim_%s' %sim)
    
    
    # ALTERNATIVE TEST
    
    # Now, we will turn off all gamma-ray components of the sky 
    # except for the component from atomic (H2) emission, as predicted 
    # from the CO12 distribution (the bad way of predicting H2) + the 
    # excess template predicted by using CO13.
    # i.e. rho(H2) = rho(CO12) * const. 

    # We can keep the default normalization of 1 to start, or we can start with 
    # some gaussian-random initial normalization.
    gta.set_norm("galdiff04",0.8) # CO13 r=0-5
    gta.set_norm("galdiff05",0.8) # CO13 r=6-9
    gta.set_norm("galdiff06",1.2) # CO13 r=10-12
    gta.set_norm("galdiff07",1.0) # CO13 r=13-16
    gta.set_norm("galdiff08",1.0) # CO13 excess template map
    gta.set_norm("galdiff09",0.0) # Smart excess template map

    
    # We could switch from the CO13 excess template to the 
    # Convolutional neural network predicted excess template
    # by turning off and fixing the modified excess template 
    # and turning on and freeing the smart exess template.

    #gta.model_counts_map("MapSource") # Write out counts map for excess template alone.
    #gta.write_model_map("ExcessCounts",name="MapSource")

    # So, free all other parameters and fit:
    gta.free_sources(free=True)
    gta.free_index("galdiff08", free=False)
    gta.free_index("galdiff09", free=False)
    gta.free_source("galdiff07",free=False)
    gta.free_source("galdiff09", free=False)

    Fit2 = gta.fit()
    alternative = Fit2["loglike"]
    gta.write_roi('after_alternative_fit_sim_%s' %sim)
    # Save the best-fit map to a file.
    gta.write_model_map("alternative_model_sim_%s" %sim)
    
    # Calculate source spectrum:
    ltcube = '/pub/abrought/fermi-ml/data/ltcube.fits'
    obs = BinnedObs(srcMaps='srcmap_00.fits',expCube=ltcube,binnedExpMap='bexpmap_00.fits',irfs='P8R3_CLEAN_V2')
    like = BinnedAnalysis(obs,'after_alternative_fit_sim_%s_00.xml' %sim, optimizer='MINUIT')
    Elist,Flist = CalcFlux(like,'galdiff08')
    data = {"energ[MeV]":Elist,"flux[MeV/cm^2/s]":Flist}
    df = pd.DataFrame(data=data)
    df.to_csv("%s/spectraout/excess_flux_sim_%s.dat" % (%indir, %sim),sep="\t",index=False)
    
    # Calculte TS:
    TS = -2*(null - alternative)

    # Write final TS to file:
    savefile = "%s/tsout/TS_sim_%s.txt" % (%indir, %sim)
    f = open(savefile,"w")
    f.write(str(TS))
    f.close()
    
    # Remove baseline ROI files to reduce storage:
    # 
    # The simulations generate a lot of memory, and the bulk of the memory is contained in the files
    # that *don't* need to get generated each time a simulation is run.
    # 
    # If you are running this for the first time, make sure to comment out this line so that the files
    # can be initially saved. Then move these files from the output/sim_0 directory to the workdir.
    # All subsequent simulations can be run with this line uncommented to clean up space after the simulation
    # is finished.
    os.system('find . -maxdepth 1 -type f -not \( -name \'*_sim_*\' -or -name \'*ccube*.fits\' -or -name \'srcmap*.fits\' \) -delete')
    
    return

################################################
if __name__=="__main__":
    main(sys.argv)
