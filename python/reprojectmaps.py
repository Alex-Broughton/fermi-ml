#Written by Chris Karwin, May 2019, UCI
#Purpose: reproject GALPROP output maps and add energy extension.
#The energy extension is needed for Fermipy.

#imports:
from astropy.io import fits
from reproject import reproject_interp
import numpy as np
import matplotlib.pyplot as plt

#the reprojection is done for each map corresponding to number of GALPROP radial bins:
for i in range(0,17):
    print(i)
    this_dir  = "/pub/abrought/fermi-ml/data/M31/"
    that_dir  = "/pub/abrought/fermi-ml/data/M31/reprojected/"
    this_file = "pion_decay_HIR_mapcube_comp_%i_56_M31_AIC.gz" %i
    savefile = that_dir + 'reprojected_' + this_file
    
    fits.setval(this_dir + this_file, 'CUNIT3', value='MeV')
    print("Fixed header...")
    
    hdu = fits.open(this_dir + this_file)[0]
    new_header = hdu.header.copy() 
    
    new_header['CRVAL1'] = 355
    new_header['CDELT1'] = -0.03125
    new_header['CRPIX1'] = 0

    print("Reprojecting...")
    new_image, footprint = reproject_interp(hdu, new_header)
    print(np.sum(new_image.ravel()))


    print("Saving...")
    fits.writeto(savefile, new_image, new_header, overwrite=True)

    #add energy extenstion:
    hdu = fits.open(this_dir + this_file)[1]
    energy = hdu.data
    energy = np.array(energy, dtype=(np.record, [('Energy', '>f8')]))
    
    hdu_2 = fits.open(savefile)
    hdu_2.append(fits.BinTableHDU(energy))
    header_1 = hdu_2[1].header
    header_1["TTYPE1"] = 'Energy'
    header_1["TUNIT1"] = 'MeV'
    header_1["EXTNAME"] = 'ENERGIES'
    hdu_2.writeto(savefile,overwrite=True)
