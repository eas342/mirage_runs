from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np

dat = ascii.read('ck_gto_f_wasp80b_1trans_R100.div',
                 names=['wavel','rprs_noise','rprs_truth','err'])

outDat = Table()
outDat['wave (um)'] = dat['wavel']
outDat['depth'] = dat['rprs_truth']
outDat['Rp/R*'] = np.sqrt(outDat['depth'])
outDat.write('ck_gto_f_wasp80b_1trans_R100.csv')