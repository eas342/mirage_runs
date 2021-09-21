from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
import sys

def reformat(filename):
    dat = ascii.read(filename,
                     names=['wavel','depth_noise','depth_truth','err'])
    if "trans_R100" in filename:
        modelType = 'transit'
    elif "eclipse_R100" in filename:
        modelType = 'eclipse'
    else:
        raise Exception("Unrecognized file type")
    
    outname = filename.replace('.div','.csv')
    outDat = Table()
    outDat['wave (um)'] = dat['wavel']
    outDat['depth'] = dat['depth_truth']
    
    if modelType == 'transit':
        outDat['Rp/R*'] = np.sqrt(outDat['depth'])
    
    outDat.write(outname)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        filename = 'ck_gto_f_wasp80b_1trans_R100.div'
    else:
        filename = sys.argv[1]
    reformat(filename)
    