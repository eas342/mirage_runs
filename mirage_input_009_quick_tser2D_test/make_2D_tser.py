import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt

def make_2D():
    
    nWaves = 10
    nTime = 100
    waves = np.linspace(2.1,4.3,nWaves)
    timeArr = np.arange(nTime)
    flux2D = np.ones([nTime,nWaves])
    
    ## make a flux bump in the middle
    flux2D[40:60,3:6] = 1.2
    
    primHDU = fits.PrimaryHDU(flux2D)
    primHDU.name = "FLUX"
    primHDU.header['AXIS1'] = ('wavelength', 'wavelength axis')
    primHDU.header['AXIS2'] = ('time','time axis')
    primHDU.header['BUNIT'] = ('norm-flux', 'normalized flux')
    timeHDU = fits.ImageHDU(timeArr)
    timeHDU.name = 'TIME'
    timeHDU.header['BUNIT'] = ('seconds','seconds')
    waveHDU = fits.ImageHDU(waves)
    waveHDU.name = 'WAVE'
    waveHDU.header['BUNIT'] = ('um','microns')
    
    
    HDUList = fits.HDUList([primHDU,timeHDU,waveHDU])
    HDUList.writeto('ex_timeser2D_longer.fits',overwrite=True)

if __name__ == "__main__":
    make_2D()
    