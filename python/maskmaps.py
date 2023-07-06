#imports:
from astropy.io import fits
from astropy.wcs import WCS
import astropy.wcs.utils as utils
import numpy as np

#upload footprint:
hdu = fits.open("/pub/abrought/fermi-ml/data/foot_print.fits")
foot_header = hdu[0].header
footprint = hdu[0].data
foot_wcs = WCS(foot_header)

#mask footprint pixels:
mask = np.argwhere(footprint > 1)

#make l and b lists: 
blist = []
llist = []
for each in mask:
    b = each[0]
    l = each[1]
    blist.append(b)
    llist.append(l)

#transfer pixels in foot_wcs to sky coordinates:
coords = utils.pixel_to_skycoord(llist,blist,foot_wcs)

#upload GALPROP:
combine = [[0,5],[6,9],[10,12],[13,16]]

for each in combine:

    low = each[0]
    high = each[1]
    
    print(low,high)
    
    this_dir  = "/pub/abrought/fermi-ml/data/M31/combined/"
    that_dir  = "/pub/abrought/fermi-ml/data/M31/masked/"
    this_file = this_dir + "reprojected_pion_decay_HIR_mapcube_comp_%s_%s_56_M31_AIC.fits" % (str(low),str(high))

    savefile = that_dir + "masked_" + "reprojected_pion_decay_HIR_mapcube_comp_%s_%s_56_M31_AIC.fits" % (str(low),str(high))

    hdu = fits.open(this_file)
    gal_header = hdu[0].header
    gal_data = hdu[0].data
    gal_wcs = WCS(gal_header)

    gal_shape = gal_data.shape

    #convert sky coordinates to pixels in gal_wcs:
    pixs =  utils.skycoord_to_pixel(coords,gal_wcs)

    mask_index = np.array(pixs)
    mask_index = mask_index.astype(int)

    #apply mask:
    for E in range(0,21):
        print("Energy bin #" + str(E))

        print("Masking footprint...")
        #footprint mask:
        for i, px in enumerate(mask_index.T):
            try:
                gal_data[E,px[1],px[0]] = 0.0
                #gal_data[E,mask_index[1],mask_index[0]] = 0.0
            except:
                continue
                
        print("Masking out of field regions...")
        #mask top and bottom region and bad field:
        for l in range(0,gal_shape[2]):
            for b in range (0,gal_shape[1]):
                this_world = gal_wcs.all_pix2world(np.array([[l,b,1]]),0)
                this_l = this_world[0][0]
                this_b = this_world[0][1]
                if this_l < 300:
                    gal_data[E,b,l] = 0
                if this_l > 350:
                    gal_data[E,b,l] = 0
                if this_b < -0.5 or this_b > 0.5:
                    gal_data[E,b,l] = 0
                if 326 <= this_l <= 327:
                    gal_data[E,b,l] = 0

    #savefile:
    hdu.writeto(savefile,overwrite=True)
