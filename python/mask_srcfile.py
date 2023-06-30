#imports:
import numpy as np
from astropy.wcs import WCS
import astropy.wcs.utils as utils
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits

hdu = fits.open("srcmap_00.fits")
header = hdu[0].header
wcs = WCS(header)

comp_list = [0,3,4,5,6,7,8,9,10,11]

mask_region = np.arange(301,351,2) # Select the alt testing region
#mask_region = np.arange(300,350,2) # Select the tiled testing region
#make mask:
for i in comp_list:
    print(i)

    data = hdu[i].data
    shape = data.shape

    for E in range(0,shape[0]):

        #mask top and bottom region and bad field:
        for l in range(0,shape[2]):
            for b in range (0,shape[1]):
                this_world = wcs.all_pix2world(np.array([[l,b,1]]),0)
                this_l = this_world[0][0]
                this_b = this_world[0][1]
                for each in mask_region:
                    low = each
                    high = low+1
                    if low < this_l < high:
                        data[E,b,l] = 0

hdu.writeto("srcmap_testing.fits",overwrite=True)
