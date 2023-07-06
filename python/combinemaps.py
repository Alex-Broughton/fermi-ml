#imports:
from astropy.io import fits

#annuli to combine:
combine = [[0,5],[6,9],[10,12],[13,16]]

for each in combine:
    
    low = each[0]
    high = each[1]
    
    this_dir  = "/pub/abrought/fermi-ml/data/M31/reprojected/"
    that_dir  = "/pub/abrought/fermi-ml/data/M31/combined/"
    this_file = this_dir + "reprojected_pion_decay_HIR_mapcube_comp_%s_56_M31_AIC.gz" % str(low)
    #this_file = "reprojected_pion_decay_H2R_mapcube_comp_%s_56_Mopra.gz" % str(low)
    
    hdu = fits.open(this_file)
    data = hdu[0].data

    for i in range(low+1,high+1):

        tmp_file = this_dir + "reprojected_pion_decay_HIR_mapcube_comp_%s_56_M31_AIC.gz" %i
        add_hdu = fits.open(tmp_file)
        add_data = add_hdu[0].data
        data += add_data
        #tmp_file.close()

    hdu.writeto(that_dir + "reprojected_pion_decay_HIR_mapcube_comp_%s_%s_56_M31_AIC.fits" % (str(low),str(high)), overwrite=True)
    #hdu.close()


