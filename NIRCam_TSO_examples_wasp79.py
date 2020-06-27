#!/usr/bin/env python
# coding: utf-8

# # Create NIRCam TSO Simulated Data

# This notebook shows how to run [Mirage](https://mirage-data-simulator.readthedocs.io/en/latest/) to create TSO data.
# 
# It uses an xml and a pointing file exported from an [APT](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview) file for a proposal containing both Grism and Imaging mode time series observations. From these files, along with several other user-inputs, Mirage will create a series of input yaml files. One yaml file is produced for each detector in each exposure. Mirage then creates the simulated data corresponding to each yaml file.
# 
# Note that in a TSO exposure, users will typically get data from 3 [detectors](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview) (A1 and A3 for shortwave data, and A5 for longwave data) in each exposure.
# 
# In this example, we will use an example APT program that looks at WASP-79. The proposal contains 3 observations:
# 
# * Grism TSO - F444W, with accompanying Weak Lens +8 and F182M shortwave imaging observation
# * Grism TSO - F322W2, with accompanying Weak Lens +8 and F210M shortwave imaging observation
# * Imaging TSO - F470N in the longwave channel, and Weak Lens +8 and F210M in the shortwave channel
# 
# We will create yaml files for all of these observations (including the accompanying Target Acquisition exposures). We will then use Mirage to create simulated data for a longwave Grism Time Series observation, one of the Imaging observations that accompany the Grism data, as well as a standalone Imaging Time Series Observation.
# 
# As written, the notebook can be executed from any directory and will find the example xml and pointing files in the Mirage repository. Outputs are all saved in the working directory.

# ## Table of Contents

# * [Imports](#imports)
# * [Inputs](#inputs)
#   * [Stellar Spectrum](#stellar_spectrum)
#   * [Batman parameters](#batman_parameters)
#   * [Transmission Spectrum - needed for GrTSO](#transmission_spectrum)
#   * [Grism TSO catalog](#grism_tso_catalog)
#   * [Lightcurve File - needed for Imaging TSO](#lightcurve_file)
#   * [Imaging TSO catalog](#imaging_tso_catalog)
#   * [Catalog of Background Sources](#background_catalog)
#   * [Create Input Yaml Files](#yaml_files)
# * [Create Simulated Data](#create_simulated_data)
#   * [Grism TSO](#grism_data)
#   * [Accompanying Imaging TSO](#accompanying_imaging_data)
#   * [Imaging TSO](#imaging_data)
# * [Run Calibration Pipeline](#calibration_pipeline)

# ## Set Environment Variables

# If you have not yet set your CRDS-related environment variables, do that here. This must be done prior to importing the CRDS software package. Your MIRAGE_DATA environment variable should also be pointing to the location of Mirage's collection of reference files.

# In[1]:


import os
#os.environ["MIRAGE_DATA"] = "/my_files/jwst/simulations/mirage_data"
os.environ["CRDS_PATH"] = os.path.expandvars("$HOME/crds_cache")
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"


# <a id="imports"></a>
# ## Imports

# In[2]:


from astropy.io import fits, ascii
from astropy.table import Table
from astropy.visualization import simple_norm, imshow_norm
from astropy import units as u
import batman
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cmx
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement
from synphot import units
import yaml

from mirage.catalogs.hdf5_catalog import save_tso
from mirage.catalogs.catalog_generator import GrismTSOCatalog, ImagingTSOCatalog, PointSourceCatalog
from mirage.catalogs.catalog_generator import TSO_GRISM_INDEX
from mirage.grism_tso_simulator import GrismTSO
from mirage.imaging_simulator import ImgSim
from mirage.seed_image.catalog_seed_image import Catalog_seed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator


# Define paths to help organize inputs and outputs

# In[3]:


path = os.path.dirname(yaml_generator.__file__)


# In[4]:


#input_data_path = os.path.abspath(os.path.join(path, '../../examples/tso_example_data'))
input_data_path = '/home/schlawin/other_code/mirage/examples/tso_example_data'


# In[5]:


input_data_path


# In[6]:


output_dir = './'
output_yaml_dir = os.path.abspath('./yaml_files')
ensure_dir_exists(output_yaml_dir)
output_data_dir = os.path.abspath('./sim_data')
ensure_dir_exists(output_data_dir)


# The yaml files that will serve as the inputs to Mirage will be saved in this directory

# In[7]:


output_yaml_dir


# The simulated data produced by Mirage will be saved in this directory

# In[8]:


output_data_dir


# Function for viewing simulated data

# In[9]:


def show(array, title, min=0, max=1000):
    """Quick view of an array.
    
    Parameters
    ----------
    array : numpy.ndimage
        2D array
        
    title : str
        Title to place on image
        
    min : int
        Signal level corresponding to bottom of color scale
        
    max : int
        Signal corresponding to top of color scale
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    norm = simple_norm(array, stretch='log', min_cut=min, max_cut=max)
    cax = ax.imshow(array, norm=norm, origin='lower')
    cbar = fig.colorbar(cax)
    plt.title(title)
    plt.show()


# <a id="prepare_inputs"></a>
# ## Prepare Inputs

# Prior to simulated data creation, there are a number of input parameters to set and files to be created. To start with, define the xml and pointing files associated with the APT file. These files can be saved from within APT by selecting them from table that appears when you select **Export** from the **File** menu. 

# In[10]:


xml_file = os.path.join(input_data_path, 'wasp-79_example_TSO.xml')
pointing_file = xml_file.replace('.xml', '.pointing')


# <a id="stellar_spectrum"></a>
# ### Stellar spectrum

# This is the spectrum of the unocculted star associated with the TSO object. For this example, we'll use the Castelli & Kurucz models in stsynphot to generate a spectrum that is similar to WASP-79. You can generate your spectrum using any tool you like, as it will eventually be saved in an hdf5 file below.

# In[11]:


t_eff = 6750  # surface temperature
metallicity = 0.03
log_g = 4.26  # surface gravity = 182 m/s^2
sp = stsyn.grid_to_spec('ck04models', t_eff, metallicity, log_g) 


# Normalize the spectrum to be k = 9.06. There are two ways you can scale your spectrum. This first is shown here, where you manually scale the spectrum before saving to the hdf5 file. The second way is to leave the scaling to Mirage. In that case, you save the spectrum as-is, and set the flux units in the hdf5 file to 'normalized'. With that option, Mirage will automatically scale the spectrum to the magnitude indicated in the grism TSO source catalog. If you choose to scale the spectrum manually (and use flux units of 'flam' in the hdf5 file), then Mirage will ignore the source magnitude listed in the grism TSO catalog, and use the saved spectrum with no changes.

# In[12]:


# Normalize the spectrum
bp = SpectralElement.from_filter('johnson_k')
vega = SourceSpectrum.from_vega()
sp_norm = sp.normalize(9.06 * units.VEGAMAG, bp, vegaspec=vega)


# Get wavelengths and flux densities of the spectrum

# In[13]:


wavelengths = sp_norm.waveset.to(u.micron)
fluxes = sp_norm(wavelengths, flux_unit='flam')


# Examine the spectrum in the wavelengths of interest

# In[14]:


f, a = plt.subplots()
a.plot(wavelengths, fluxes)
a.set_xlim(1, 5)
a.set_ylim(0, 1e-13)
a.set_xlabel('Wavelength (microns)')
a.set_ylabel('Flux density (FLAM)')
plt.show()


# Set the units for the wavelength and flux density arrays. It's generally recommended to use flux denisty units of FLAM (erg / s / cm^2 / 𝐴˚). 

# In[15]:


wavelength_units = 'microns'
flux_units = 'flam'
#flux_units = 'normalized'


# Name of the file in which to save the spectrum

# In[16]:


sed_file = os.path.join(output_dir, 'test_grism_tso_sed_file_wasp79.hdf5')


# HDF5 files can contain multiple datasets, so make the flux and wavelength variables into lists, even though for TSO observations there will be only one source.

# In[17]:


fluxes = [fluxes]
wavelengths = [wavelengths]


# Save the spectrum in the hdf5 file. Give the dataset a special index number, defined by TSO_GRISM_INDEX, to help Mirage keep this source separate from any provided by other catalogs.

# In[18]:


with h5py.File(sed_file, "w") as file_obj:
    for i in range(len(fluxes)):
        dset = file_obj.create_dataset(str(i+TSO_GRISM_INDEX), data=[wavelengths[i].value, fluxes[i].value],
                                       dtype='f', compression="gzip", compression_opts=9)
        dset.attrs[u'wavelength_units'] = wavelength_units
        dset.attrs[u'flux_units'] = flux_units


# <a id="batman_parameters"></a>
# ### Batman Parameters

# Model lightcurves are generated using the [Batman](https://github.com/lkreidberg/batman) package. Here, set the parameters that Batman needs. See the [Batman documentation](https://www.cfa.harvard.edu/~lkreidberg/batman/) for details.

# In[19]:


# In this case, we'll use parameters based on WASP-79b, but modified for a shorter exposure time
# and a deeper eclipse, in order to save processing time and make the resulting data easier to
# examine.
params = batman.TransitParams()       # object to store transit parameters 
params.t0 = 280.                        # time of inferior conjunction 
params.per = 3162.24                       # orbital period 
params.rp = 0.723                       # planet radius (in units of stellar radii) 
params.a = 9.37                        # semi-major axis (in units of stellar radii) 
params.inc = 83.3                      # orbital inclination (in degrees) 
params.ecc = 0.                       # eccentricity 
params.w = 90.                        # longitude of periastron (in degrees) 
params.limb_dark = "nonlinear"        # limb darkening model 
params.u = [0.5, 0.1, 0.1, -0.1]      # limb darkening coefficients [u1, u2, u3, u4] 


# Generate the array of times at which the lightcurve will be calculated. In this case we know that the total exposure time is about 570 seconds, so we extend the array of times just beyond that range.

# In[20]:


times = np.linspace(0, 580, 1000)  # times at which to calculate light curve 


# <a id="transmission_spectrum"></a>
# ### Create a Transmission Spectrum

# The transmission spectrum is the wavelength-dependent effective radius of the planet, in units of the stellar radius. It must be saved in an ascii file, and will be used in the creation of the Grism TSO data. For this example, we use the simple case of a flat spectrum. 

# In[21]:


waves = np.linspace(0.9, 5.5, 1000)  # microns
trans = np.repeat(params.rp, 1000)  # R_planet / R_star


# In[22]:


tran_spec_file = os.path.join(output_dir,'transmission_spectrum.txt')
tab = Table()
tab['Wavelength'] = waves
tab['Transmission'] = trans
tab.write(tran_spec_file, format='ascii', overwrite=True)


# Plot transmission spectrum

# In[23]:


f, a = plt.subplots()
a.scatter(waves, trans, color='red', marker='o')
a.set_xlabel('Wavelength (microns)')
a.set_ylabel('Transmission')
plt.show()


# <a id="grism_tso_catalog"></a>
# ### Create Grism TSO catalog

# With the stellar spectrum and transmission spectrum saved, we can now create the Grism TSO source catalog that Mirage will use when creating the simulated data.

# In[24]:


# Name of catalog file to hold information on the TSO source
grism_tso_catalog = os.path.join(output_dir,'tso_grism_source.cat')


# Basic information on the source. Note that the object magnitude will be ignored if the saved stellar spectrum is in units of FLAM. Conversely, if the stellar spectrum has units of 'normalized', then Mirage will scale the spectrum to the magnitude indicated below.

# In[25]:


object_ra = 66.37090333
object_dec = -30.60044722
object_f444w_mag = 9.0
object_f322w2_mag = 9.05


# Create the Grism TSO catalog object and populate RA, Dec, Batman parameters, times, and the name of the transmission spectrum file.

# In[26]:


grism_cat = GrismTSOCatalog(ra=[object_ra], dec=[object_dec], semimajor_axis=[params.a],
                            orbital_inclination=[params.inc], eccentricity=[params.ecc],
                            orbital_period=[params.per], longitude_of_periastron=[params.w],
                            limb_dark_model=[params.limb_dark], limb_dark_coeffs=[params.u],
                            time_units=['second'], start_time=[np.min(times)],
                            end_time=[np.max(times)], inferior_conj=[params.t0],
                            transmission_spectrum=[tran_spec_file])


# Add the source magnitudes to the catalog

# In[27]:


# Add source magnitude
grism_cat.add_magnitude_column([object_f444w_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f444w')
grism_cat.add_magnitude_column([object_f322w2_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f322w2')


# In[28]:


grism_cat.save(grism_tso_catalog)


# Examine the contents of the catalog. There should be only a single source.

# In[29]:


grism_cat.table


# <a id="lightcurve_file"></a>
# ### Create Lightcurve File

# For the imaging time series observations, we need to provide a file that contains the lightcurve to use when simulating the data.

# Initialize the batman model using the parameters specified in the [Batman Parameters](batman_parameters) section and generate a lightcurve.

# In[30]:


m = batman.TransitModel(params, times)
flux = m.light_curve(params)


# Plot the lightcurve to be used to generate the data

# In[31]:


f, a = plt.subplots()
a.scatter(times, flux, color='red', marker='v')
a.set_xlabel('Time (sec)')
a.set_ylabel('Normalized Signal')
plt.show()


# In[32]:


lightcurve_file = os.path.join(output_dir, 'example_lightcurve.hdf5')


# Place the lightcurve into a dictionary to prepare for saving. The keys are object indexes corresponding to objects in the Mirage input catalogs. 

# In[33]:


contents = {}
contents['1'] = {'times': times,
                 'fluxes': flux}


# In[34]:


# Save
save_tso(contents, lightcurve_file, time_unit='second')


# <a id="imaging_tso_catalog"></a>
# ### Create Imaging TSO catalog

# In this case, rather than the `tso_grism_catalog` in the yaml file, the user must supply a `tso_imaging_catalog`. This catalog will contain the list of sources whose flux will be varying with time. As in the grism TSO case, objects listed in the other types of catalogs will be added to the simulation in order to create a more realistic scene.
# 
# See the [notebook on catalog creation](https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb) in the examples directory of the Mirage repository for examples of how to create source catalogs. Catalogs may also be created manually.
# 
# Mirage does not yet support the use of input spectra in order to determine filter-based magnitudes for imaging simulations. Until that ability is added, users must estimate the source magnitude in each filter.

# In[35]:


imaging_tso_catalog = os.path.join(output_dir, 'tso_imaging_source.cat')


# In[36]:


tsimg_cat = ImagingTSOCatalog(ra=[object_ra], dec=[object_dec], lightcurve_file=[lightcurve_file])


# In[37]:


object_f182m_mag = 10.0
object_f210m_mag = 9.5
object_f470n_mag = 9.0 


# In[38]:


# Add source magnitudes
tsimg_cat.add_magnitude_column([object_f182m_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f182m')
tsimg_cat.add_magnitude_column([object_f210m_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f210m')
tsimg_cat.add_magnitude_column([object_f470n_mag], magnitude_system='vegamag',
                               instrument='nircam', filter_name='f470n')


# In[39]:


tsimg_cat.save(imaging_tso_catalog)


# In[40]:


tsimg_cat.table


# <a id="background_catalog"></a>
# ### Catalog of Background Sources

# Create a catalog containing nearby sources. Mirage will add these to the data. This can help inform desired roll angles in order to avoid contamination.
# 
# Mirage accepts a series of ascii source catalogs. In this case, the options that can be used are:
# 
# 1. Point source catalog
# 2. Galaxy catalog
# 3. Extended source catalog
# 
# The point source, galaxy, and extended source catalogs contain "background" sources. That is, sources other than the primary TSO target. All catalogs are input to Mirage as entries within the input yaml file.
# 
# See the notebook on [catalog creation](https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb) in the examples directory of the Mirage repository for examples of how to create source catalogs. Catalogs may also be created manually. For this example, we will create a point source catalog with a single star.

# In[41]:


bkgd_sources_ra = [66.378244]
bkgd_sources_dec = [-30.598312]
bkgd_sources_f444w_mag = [12.50]
bkgd_sources_f322w2_mag = [12.75]
bkgd_sources_f182m_mag = [13.75]
bkgd_sources_f210m_mag = [13.39]
bkgd_sources_f470n_mag = [12.40]


# In[42]:


bkgd_cat_file = os.path.join(output_dir, 'ptsrcs.cat')


# In[43]:


bkgd_cat = PointSourceCatalog(ra=bkgd_sources_ra, dec=bkgd_sources_dec)


# In[44]:


# Add source magnitudes
bkgd_cat.add_magnitude_column(bkgd_sources_f182m_mag, magnitude_system='vegamag',
                               instrument='nircam', filter_name='f182m')
bkgd_cat.add_magnitude_column(bkgd_sources_f444w_mag, magnitude_system='vegamag',
                               instrument='nircam', filter_name='f444w')
bkgd_cat.add_magnitude_column(bkgd_sources_f322w2_mag, magnitude_system='vegamag',
                               instrument='nircam', filter_name='f322w2')
bkgd_cat.add_magnitude_column(bkgd_sources_f210m_mag, magnitude_system='vegamag',
                               instrument='nircam', filter_name='f210m')
bkgd_cat.add_magnitude_column(bkgd_sources_f470n_mag, magnitude_system='vegamag',
                               instrument='nircam', filter_name='f470n')


# In[45]:


bkgd_cat.save(bkgd_cat_file)


# In[46]:


bkgd_cat.table


# <a id="yaml_files"></a>
# ### Create input yaml files for Mirage

# Mirage requires one input yaml file for each detector and exposure. Using the files exported from ATP combined with the various files produced above, Mirage can automatically generate all of the yaml files necessary to produce all of the data in the APT file. See the Mirage documentation on [creating yaml files from APT inputs](https://mirage-data-simulator.readthedocs.io/en/latest/yaml_generator.html) for details on the input parameters used below.

# Populate the input catalog dictionary with the names of the source catalogs created above.

# In[47]:


catalogs = {'WASP-79': {'point_source': bkgd_cat_file,
                        'tso_imaging_catalog': imaging_tso_catalog,
                        'tso_grism_catalog': grism_tso_catalog,
                        }
           }


# Set the desired background level for the observations. Most common will to specify 'low', 'medium', or 'high'. The definitions of these three levels match those used in the ETC. Note that for the grism observations, if the `dateobs_for_background` is set to True, the background will be created based on the given observation date rather than the level given below.

# In[48]:


background = 'medium'


# Set the desired telescope roll angle and date for the observation. 

# In[49]:


pav3 = 0.
dates = '2021-10-25'


# In[50]:


yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                              output_dir=output_yaml_dir, simdata_output_dir=output_data_dir,
                              background=background, roll_angle=pav3,
                              dates=dates, datatype='linear, raw', dateobs_for_background=True,
                              reffile_defaults='crds')

yam.use_linearized_darks = True
yam.create_inputs()


# List the yaml files created by the generator

# In[51]:


yam.yaml_files


# Create a table showing details of what exposure each yaml file describes

# In[52]:


basenames = []
modes = []
apertures = []
filters = []
pupils = []
for yfile in yam.yaml_files:
    with open(yfile) as file_obj:
        info = yaml.safe_load(file_obj)
    basenames.append(os.path.basename(yfile))
    modes.append(info['Inst']['mode'])
    apertures.append(info['Readout']['array_name'])
    filters.append(info['Readout']['filter'])
    pupils.append(info['Readout']['pupil'])
info_tab = Table()
info_tab['Filename'] = basenames
info_tab['Mode'] = modes
info_tab['Aperture'] = apertures
info_tab['Filter'] = filters
info_tab['Pupil'] = pupils


# Exposures using `imaging` mode are for the target acquisition exposures. Those with `ts_grism` mode describe the Grism time series observations. Those with `ts_imaging` mode are for Imaging time series observations. This includes the data collected with the shortwave channel detectors while the longwave detector is producing the grism time series observations.

# In[53]:


info_tab


# <a id="create_simulated_data"></a>
# ## Create Simulated Data

# Simulate data from some of the files above.

# <a id="grism_data"></a>
# ### Grism TSO Data

# First let's simulate a grism time series exposure using the F444W filter. 

# In[53]:


gr_tso_yaml_file = os.path.join(output_yaml_dir, 'jw88888001001_01101_00002_nrca5.yaml')


# In[54]:


gr_tso_yaml_file


# In[55]:
if os.path.exists('WFSS_background.fits'):
    os.replace('WFSS_background.fits','WFSS_background_old.fits')

# NOTE: This cell will take a while (~30 min) to run 
gr_f444w = GrismTSO(gr_tso_yaml_file, SED_file=sed_file, SED_normalizing_catalog_column=None,
                    final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                    extrapolate_SED=True, override_dark=None, disp_seed_filename=None,
                    orders=["+1", "+2"])
gr_f444w.create()


# #### Examine the results

# Mirage will split the exposure into mulitple files if it is too large, in the same way that the JWST calibration pipeline will do. Exposures will be split around integrations. That is, an integration will not be split between two files.

# First, look at the seed image, which is the noiseless scene that contains only signal from the astronomical sources, including background.

# In[56]:


gr_f444w.seed_files


# Read in the first seed file. Checking the array dimensions, we see that it holds the first 60 integrations.

# In[57]:


data_444 = fits.getdata(gr_f444w.seed_files[0])


# In[58]:


data_444.shape


# Look at the seed image. The TSO source is the lower source of the two. The upper source is the background star contained in the point source catalog.

# In[59]:


show(data_444[0, 4, :, :], 'F444W: Final group of integration 0')


# Using the seed image, plot the signal along the trace of the TSO object for each integration, as a sanity check to see that the signal does change with time.

# In[60]:


def multi_spec_plot(array, yval, xstart):
    """Plot the signal along the `yval` row from `xstart` to the end,
    in the final group of each integration. This will be a quick
    sanity check to show how the signal is changing with wavelength
    across the integrations. Step along the color table and use a
    different color for each integration.
    
    Parameters
    ----------
    array : numpy.ndarray
        4D array containing the spectrum in each group and integration
        
    yval : int
        Row number to plot
        
    xstart : int
        First column number to plot. Plot from here to the end of the row
    """
    num_integrations = array.shape[0]
    num_grps = array.shape[1]
    clrs = cmx.rainbow(np.linspace(0, 1, num_integrations))
    f, a = plt.subplots()
    f.set_size_inches(11.5,9)
    for integration, color in zip(range(num_integrations), clrs): 
        slopes = array[integration, num_grps - 1, yval, xstart:] 
        a.plot(np.arange(len(slopes)), slopes, color=color)
        a.set_xlabel('Column Number')
        a.set_ylabel('Signal (DN)')


# The plotting function above uses the `rainbow` color table, which starts with blue and moves through the spectrum to red. Below we plot the signal in the 51st row of the aperture in the final group of each of the 60 integrations. We see that the initial integrations (blue, but covered and not visible) and the final integrations (red) overlay each other, as expected since they occur before and after the transit, respectively. In the middle of the exposure, the signal in the spectrum decreases with time (light blue -> green) down to a minimum and then increases again (orange), as expected. The F444W filter covers a wavelength range of approximately 3.8 - 5.1 microns.

# In[61]:


multi_spec_plot(data_444, 51, 700)


# Look at the raw, or uncalibrated data

# In[62]:


uncal_file = os.path.join(output_data_dir, 'jw88888001001_01101_00002-seg001_nrca5_uncal.fits')


# In[63]:


uncal = fits.getdata(uncal_file)


# Look at a difference image between the last and first groups of an integration. By using the difference, we remove bias offsets and make the sources more visible. Raw files have signals in integers, so multiply by 1.0 in order to make them floats before subtracting. In this case, the signal from the background star is faint enough that it is not visible here above the noise.
# 
# The elevated background signal in the 512-column wide area is the result of each 512-column block being read out by a different amplifier. These differences will be minimized once the data have been run through the calibration pipeline.

# In[64]:


diff = 1. * uncal[0, 4, :, :] - 1. * uncal[0, 0, :, :]


# In[65]:


show(diff, "Uncalibrated Data", min=0, max=400)


# Create simulated data for the same observation, but in the F322W2 filter

# In[54]:


gr_tso_f322w2_yaml_file = os.path.join(output_yaml_dir, 'jw88888002001_01101_00002_nrca5.yaml')


# Rename the background file
os.replace('WFSS_background.fits','WFSS_background_F444W.fits')

# In[ ]:


gr_f322w2 = GrismTSO(gr_tso_f322w2_yaml_file, SED_file=sed_file, SED_normalizing_catalog_column=None,
                     final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                     extrapolate_SED=True, override_dark=None, disp_seed_filename=None,
                     orders=["+1", "+2"])
gr_f322w2.create()


# Read in one of the seed image files and examine the data



data_322 = fits.getdata(gr_f322w2.seed_files[0])


# Note how with F322W2 as the crossing filter, the location of the dispersed objects have shifted to the left relative to what they were with the F444W filter.

# In[69]:


show(data_322[0, 4, :, :], 'F322W2: Final group of integration 0')


# Plot the signal along the trace, as with the F444W data. Again we see that the initial (blue) and final (red) integrations, which were outside of the transit, show identical signals at the maximum level, while those integrations that occured within the transit (light blue, green, orange) show reduced signal. The F322w2 filter covers a wavelength range of approximately 2.4 - 4.0 microns.

# In[70]:


multi_spec_plot(data_322, 51, 50)


# <a id="accompanying_imaging_data"></a>
# ### Accompanying Imaging TSO Data

# Let's look at the imaging time series data collected on one of the shortwave detectors while the grism data above was being collected in the LW channel.

# In[71]:


img_tso_sw_yaml = os.path.join(output_yaml_dir, 'jw88888001001_01101_00001_nrca1.yaml')


# In[ ]:


img_tso = ImgSim()
img_tso.paramfile = img_tso_sw_yaml
img_tso.create()


# Examine the raw output file. This is the simulated data before any calibration. Detector A1 covers the same area as the [left half of detector A5](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview) (which was used to collect the Grism TSO above). This is why the source appears towards the right edge of A1. The +8 weak lens was also used in this observation, resulting in the expanded PSF.

# In[ ]:


img_seed_file = os.path.join(output_data_dir,
                             'jw88888001001_01101_00001_nrca1_uncal_F182M_seg001_part001_seed_image.fits')


# In[ ]:


img_seed = fits.getdata(img_seed_file)


# In[ ]:


show(img_seed[0, 4, :, :], 'Noiseless Seed Image', min=0, max=3000)


# Again, let's look at signal levels in all of the integrations, using the noiseless seed image. Since there is no noise and no background, we'll do a simple sum of the signal in a small box. Here we see that the depth of the transit matches that from the lightcurve we created with the Batman parameters above.

# In[ ]:


signals = []
cx, cy = 1857, 151
for i in range(img_seed.shape[0]):
    signals.append(np.sum(img_seed[i, 4, cy-1: cy+2, cx-1: cx+2]))


# In[ ]:


f, a = plt.subplots()
a.plot(np.arange(len(signals)), signals, color='red')
a.set_xlabel('Integration Number')
a.set_ylabel('Aperture Signal')


# Now look at the raw data

# In[ ]:


uncal_file = os.path.join(output_data_dir, 'jw88888001001_01101_00001-seg001_nrca1_uncal.fits')


# In[ ]:


data = fits.getdata(uncal_file)


# In[ ]:


data.shape


# By displaying the difference between the last and first groups of one integration, we can remove bias signal that would otherwise make the object difficult to see.

# In[ ]:


cds = 1.* data[0, 4, :, :] - 1. * data[0, 0, :, :]


# The image appears to have 4 512-column wide blocks because each block has been read out by a different amplifier. These differences will be minimized once the data have been run through the calibration pipeline. Note that this observation was taken with one of NIRCam's weak lenses in place, which is why the star's PSF has become so large.

# In[ ]:


show(cds, 'Difference Image', min=0, max=1000)


# In[ ]:





# <a id="imaging_data"></a>
# ### Imaging TSO Observation

# The imaging time series observation in this proposal contains an observation of WASP-79 using the SUB400P subarray. This is a 400x400 pixel subarray located in the upper right corner of the B1 and B5 detectors. Run Mirage on both the B1 and B5 exposures.

# In[ ]:


img_tso_sw_yaml = os.path.join(output_yaml_dir, 'jw88888003001_01101_00001_nrcb1.yaml')


# In[ ]:


img_tso = ImgSim()
img_tso.paramfile = img_tso_sw_yaml
img_tso.create()


# In[ ]:


img_tso_lw_yaml = os.path.join(output_yaml_dir, 'jw88888003001_01101_00002_nrcb5.yaml')


# In[ ]:


img_tso = ImgSim()
img_tso.paramfile = img_tso_lw_yaml
img_tso.create()


# ### Look at the seed image

# In[ ]:


seed_image_file = os.path.join(output_data_dir, 'jw88888003001_01101_00001_nrcb1_uncal_F210M_seed_image.fits')


# In[ ]:


with fits.open(seed_image_file) as hdulist:
    sw_img_data = hdulist[1].data
    grptime = hdulist[0].header['TGROUP']
sw_img = sw_img_data[0, -1, :, :]


# Look at the final group of the initial integration on the B1 detector. The edges of the PSF stamp from the library are visible against the zero-signal background.

# In[ ]:


show(sw_img,'Imaging TSO', min=0, max=1000)


# And the same image for the B5 detector

# In[ ]:


lw_seed_image_file = os.path.join(output_data_dir, 'jw88888003001_01101_00002_nrcb5_uncal_F470N_seed_image.fits')


# In[ ]:


with fits.open(lw_seed_image_file) as hdulist:
    lw_img_data = hdulist[1].data
lw_img = lw_img_data[0, -1, :, :]


# The pixel scale in the longwave (B5) detector is is twice as large as in the shortwave detector, so the same 400x400 pixel aperture covers 4 times the area on the sky compared to the B1 case above.

# In[ ]:


show(lw_img,'Imaging TSO', min=0, max=1000)


# ### Rough photometry and lightcurve

# As was done in the [Imaging TSO](#imaging_data) section, we perform rough photometry on each of the integrations within the exposure in an effort to plot the lightcurve. As before, the depth of the measured lightcurve matches that in the initial Batman parameters defined above.

# In[ ]:


num_ints, num_groups, ny, nx = sw_img_data.shape
sw_img_phot = np.zeros(num_ints)
for integ in range(num_ints):
    sw_img_phot[integ] = np.sum(sw_img_data[integ, -1, 100:300, 100:300])


# In[ ]:


f, a = plt.subplots()
a.plot(np.arange(len(sw_img_phot)) * grptime * (num_groups + 1), sw_img_phot, color='red')
a.set_xlabel('Time (sec)')
a.set_ylabel('Total Aperture Signal')


# <a id="calibration_pipeline"></a>
# ## Run the Calibration Pipeline

# In[ ]:





# See the accompanying notebook for a demonstration of how to run the calibration pipeline on these data.
