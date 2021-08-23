#!/usr/bin/env python
# coding: utf-8

# # Create NIRCam TSO Simulated Data with more general inputs

# This notebook shows how to run [Mirage](https://mirage-data-simulator.readthedocs.io/en/latest/) to create TSO data.
# 
# It uses an xml and a pointing file exported from an [APT](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview) file for a proposal containing both Grism and Imaging mode time series observations. From these files, along with several other user-inputs, Mirage will create a series of input yaml files. One yaml file is produced for each detector in each exposure. Mirage then creates the simulated data corresponding to each yaml file.
# 
# Note that in a TSO exposure, users will typically get data from 3 [detectors](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview) (A1 and A3 for shortwave data, and A5 for longwave data) in each exposure.
# 
# In this example, we will use an example APT program that looks at WASP-43. The proposal contains 3 observations:
# 
# 
# * Grism TSO - F322W2, with accompanying Weak Lens +8 and F210M shortwave imaging observation
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
#os.environ["CRDS_PATH"] = "/fenrirdata1/es_tso/crds_cache"
#os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ['MIRAGE_DATA']


# In[2]:


os.environ["CRDS_PATH"], os.environ["CRDS_SERVER_URL"]


# <a id="imports"></a>
# ## Imports

# In[3]:


from astropy.io import fits, ascii
from astropy.table import Table
from astropy.visualization import simple_norm, imshow_norm
from astropy import units as u
from astropy.io import fits, ascii
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
from mirage.catalogs import create_catalog
from mirage.grism_tso_simulator import GrismTSO
from mirage.imaging_simulator import ImgSim
from mirage.seed_image.catalog_seed_image import Catalog_seed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator


# Define paths to help organize inputs and outputs

# In[4]:


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
    plt.savefig(os.path.join(output_dir,'{}.pdf'.format(title.replace(' ','_'))))


class grismWrapper(object):
# In[5]:

    def __init__(self,input_data_path = 'mirage_input_001/',
                 output_data_dir = 'mirage_output_001',
                 xml_file="example_input.xml",
                 source_params="source_params.yaml"):
        """
        A wrapper for running NIRCam grism TSO simulations w/ MIRAGE
        
        input_data_path: str
            Directory for input data: XML, source parameters. This
            is also where Mirage will dump yaml files
        output_data_dir: str
            Directory for output data
        xml_file: str
            Name of the XML file made by APT. Put that and .pointings in the input path
        source_params: str
            The parameters of the system with the exoplanet
        """
        # In[7]:
        self.input_data_path = input_data_path
        self.output_data_dir = os.path.abspath(output_data_dir)
        print("output data dir: {}".format(self.output_data_dir))

        # The yaml files that will serve as the inputs to Mirage will be saved in this directory
        self.output_yaml_dir = os.path.abspath(os.path.join(input_data_path,'yaml_files'))
        print("output yaml: {}".format(self.output_yaml_dir))
        
        dirList = [self.input_data_path,self.output_data_dir,self.output_yaml_dir]
        for oneDir in dirList:
            ensure_dir_exists(oneDir)
        
        self.xml_file = os.path.join(input_data_path, xml_file)
        
        self.pointing_file = xml_file.replace('.xml', '.pointing')
        
        
        # <a id="stellar_spectrum"></a>
        # ### Stellar spectrum
        
        # This is the spectrum of the unocculted star associated with the TSO object. For this example, we'll use the Castelli & Kurucz models in stsynphot to generate a spectrum that is similar to WASP-79. You can generate your spectrum using any tool you like, as it will eventually be saved in an hdf5 file below.
        
        # In[13]:
        
        self.sys_param_path = os.path.join(self.input_data_path,source_params)
        with open(self.sys_param_path) as sys_f:
            self.sys_params = yaml.safe_load(sys_f)
    
    
    
    def prep_sources(self):
        
        t_eff = self.sys_params['stellar']['teff']  # surface temperature
        metallicity = self.sys_params['stellar']['metallicity']
        log_g = self.sys_params['stellar']['logg'] # surface gravity
        sp = stsyn.grid_to_spec('ck04models', t_eff, metallicity, log_g) 
        
        
        # Normalize the spectrum to be k = kmag. There are two ways you can scale your spectrum. This first is shown here, where you manually scale the spectrum before saving to the hdf5 file. The second way is to leave the scaling to Mirage. In that case, you save the spectrum as-is, and set the flux units in the hdf5 file to 'normalized'. With that option, Mirage will automatically scale the spectrum to the magnitude indicated in the grism TSO source catalog. If you choose to scale the spectrum manually (and use flux units of 'flam' in the hdf5 file), then Mirage will ignore the source magnitude listed in the grism TSO catalog, and use the saved spectrum with no changes.
        
        # In[49]:
        
        # Normalize the spectrum
        bp = SpectralElement.from_filter('johnson_j')
        vega = SourceSpectrum.from_vega()
        sp_norm = sp.normalize(self.sys_params['system']['kmag'] * units.VEGAMAG, bp, vegaspec=vega)
        
        
        # Get wavelengths and flux densities of the spectrum
        
        # In[50]:
        
        
        wavelengths = sp_norm.waveset.to(u.micron)
        fluxes = sp_norm(wavelengths, flux_unit='flam')
        
        
        # Examine the spectrum in the wavelengths of interest
        
        # In[56]:
        
        
        f, a = plt.subplots()
        a.plot(wavelengths, fluxes)
        a.set_xlim(1, 5)
        a.set_ylim(0, 2.5e-13)
        a.set_xlabel('Wavelength (microns)')
        a.set_ylabel('Flux density (FLAM)')
        f.savefig(os.path.join(output_dir,'star_spec_check.png'),dpi=150)
        
        # Set the units for the wavelength and flux density arrays. It's generally recommended to use flux denisty units of FLAM (erg / s / cm^2 / 𝐴˚). 
        
        # In[57]:
        
        
        wavelength_units = 'microns'
        flux_units = 'flam'
        #flux_units = 'normalized'
        
        
        # Name of the file in which to save the spectrum
        
        # In[58]:
        
        
        sed_file = os.path.join(output_dir, 'test_grism_tso_sed_file_wasp43.hdf5')
        
        
        # **ES added** Trim the files to just the relevant wavelengths
        
        # Show how many wavelength points were kept versus before
        
        # In[62]:
        
        
        pts_keep = (wavelengths > 2.0 * u.micron) & (wavelengths < 5.4 * u.micron)
        np.sum(pts_keep), wavelengths.shape
        
        
        # HDF5 files can contain multiple datasets, so make the flux and wavelength variables into lists, even though for TSO observations there will be only one source.
        
        # In[63]:
        
        
        output_fluxes = [fluxes[pts_keep]]
        output_wavelengths = [wavelengths[pts_keep]]
        
        
        # Save the spectrum in the hdf5 file. Give the dataset a special index number, defined by TSO_GRISM_INDEX, to help Mirage keep this source separate from any provided by other catalogs.
        
        # In[64]:
        
        
        with h5py.File(sed_file, "w") as file_obj:
            for i in range(len(output_fluxes)):
                oneDat = [output_wavelengths[i].value, output_fluxes[i].value]
                dset = file_obj.create_dataset(str(i+TSO_GRISM_INDEX), data=oneDat,
                                               dtype='f', compression="gzip", compression_opts=9)
                dset.attrs[u'wavelength_units'] = wavelength_units
                dset.attrs[u'flux_units'] = flux_units
        
        
        # <a id="batman_parameters"></a>
        # ### Batman Parameters
        
        # Model lightcurves are generated using the [Batman](https://github.com/lkreidberg/batman) package. Here, set the parameters that Batman needs. See the [Batman documentation](https://www.cfa.harvard.edu/~lkreidberg/batman/) for details.
        
        # In[65]:
        
        
        # Esposito et al. 2017 parameters on NASA Exoplanet Archive
        params = batman.TransitParams()       # object to store transit parameters 
        params.t0 = self.sys_params['planet']['t0']                        # time of inferior conjunction (This one is just offset from start)
        params.per = self.sys_params['planet']['per']    # orbital period, Esposito
        params.rp = self.sys_params['planet']['rp']                      # planet radius (in units of stellar radii) , Esposito
        params.a = self.sys_params['planet']['a']                         # semi-major axis (in units of stellar radii), Esposito
        params.inc = self.sys_params['planet']['inc']                     # orbital inclination (in degrees) , Esposito
        params.ecc = self.sys_params['planet']['ecc']                      # eccentricity 
        params.w = self.sys_params['planet']['w']                        # longitude of periastron (in degrees) 
        params.limb_dark = self.sys_params['planet']['limb_model']        # limb darkening model 
        params.u = self.sys_params['planet']['ld_u']      # limb darkening coefficients [u1, u2, u3, u4] 
        
        
        # In[66]:
        
        
        print("Transit center param={}".format(params.t0))
        
        
        # Generate the array of times at which the lightcurve will be calculated. In this case we know that the total exposure time is about 570 seconds, so we extend the array of times just beyond that range.
        
        # In[67]:
        
        
        times = np.linspace(0, self.sys_params['obs']['dur'], 1000)  # times at which to calculate light curve 
        
        
        # <a id="transmission_spectrum"></a>
        # ### Create a Transmission Spectrum
        
        # The transmission spectrum is the wavelength-dependent effective radius of the planet, in units of the stellar radius. It must be saved in an ascii file, and will be used in the creation of the Grism TSO data. For this example, we use the simple case of a flat spectrum. 
        
        # In[68]:
        
        
        dat = ascii.read(self.sys_params['planet']['spec_file'])
        dat.colnames
        
        
        # In[69]:
        
        
        waves = dat['wave (um)']
        trans = dat['Rp/R*']
        
        
        # In[70]:
        
        
        # waves = np.linspace(0.9, 5.5, 1000)  # microns
        # trans = np.repeat(params.rp, 1000)  # R_planet / R_star
        
        
        # In[71]:
        
        
        tran_spec_file = os.path.join(self.output_dir,'transmission_spectrum.txt')
        tab = Table()
        tab['Wavelength'] = waves
        tab['Transmission'] = trans
        tab.write(tran_spec_file, format='ascii', overwrite=True)
        
        
        # Plot transmission spectrum
        
        # In[72]:
        
        
        f, a = plt.subplots()
        a.scatter(waves, trans, color='red', marker='o')
        a.set_xlabel('Wavelength (microns)')
        a.set_ylabel('Transmission')
        f.savefig(os.path.join(output_dir,'transmission_spec_check.png'),dpi=150)
        
        
        
        # <a id="grism_tso_catalog"></a>
        # ### Create Grism TSO catalog
        
        # With the stellar spectrum and transmission spectrum saved, we can now create the Grism TSO source catalog that Mirage will use when creating the simulated data.
        
        # In[73]:
        
        
        # Name of catalog file to hold information on the TSO source
        grism_tso_catalog = os.path.join(output_dir,'tso_grism_source.cat')
        
        
        # Basic information on the source. Note that the object magnitude will be ignored if the saved stellar spectrum is in units of FLAM. Conversely, if the stellar spectrum has units of 'normalized', then Mirage will scale the spectrum to the magnitude indicated below.
        
        # In[74]:
        
        
        object_ra = self.sys_params['system']['ra']
        object_dec = self.sys_params['system']['dec']
        object_f444w_mag = self.sys_params['system']['f444w_mag']
        object_f322w2_mag = self.sys_params['system']['f322w2_mag']
        
        
        # Create the Grism TSO catalog object and populate RA, Dec, Batman parameters, times, and the name of the transmission spectrum file.
        
        # In[75]:
        
        
        grism_cat = GrismTSOCatalog(ra=[object_ra], dec=[object_dec], semimajor_axis=[params.a],
                                    orbital_inclination=[params.inc], eccentricity=[params.ecc],
                                    orbital_period=[params.per], longitude_of_periastron=[params.w],
                                    limb_dark_model=[params.limb_dark], limb_dark_coeffs=[params.u],
                                    time_units=['second'], start_time=[np.min(times)],
                                    end_time=[np.max(times)], inferior_conj=[params.t0],
                                    transmission_spectrum=[tran_spec_file])
        
        
        # Add the source magnitudes to the catalog
        
        # In[76]:
        
        
        # Add source magnitude
        grism_cat.add_magnitude_column([object_f444w_mag], magnitude_system='vegamag',
                                       instrument='nircam', filter_name='f444w')
        grism_cat.add_magnitude_column([object_f322w2_mag], magnitude_system='vegamag',
                                       instrument='nircam', filter_name='f322w2')
        
        
        # In[77]:
        
        
        grism_cat.save(grism_tso_catalog)
        self.grism_tso_catalog = grism_tso_catalog

    
        # Examine the contents of the catalog. There should be only a single source.
        
        # In[78]:
        
        
        grism_cat.table
        
        
        # <a id="lightcurve_file"></a>
        # ### Create Lightcurve File
        
        # For the imaging time series observations, we need to provide a file that contains the lightcurve to use when simulating the data.
        
        # Initialize the batman model using the parameters specified in the [Batman Parameters](batman_parameters) section and generate a lightcurve.
        
        # In[79]:
        
        
        m = batman.TransitModel(params, times)
        flux = m.light_curve(params)
        
        
        # Plot the lightcurve to be used to generate the data
        
        # In[80]:
        
        
        f, a = plt.subplots()
        a.scatter(times, flux, color='red', marker='v')
        a.set_xlabel('Time (sec)')
        a.set_ylabel('Normalized Signal')
        f.savefig(os.path.join(output_dir,'tseries_check.png'),dpi=150)
        
        
        # In[81]:
        
        
        lightcurve_file = os.path.join(output_dir, 'example_lightcurve.hdf5')
        
        
        # Place the lightcurve into a dictionary to prepare for saving. The keys are object indexes corresponding to objects in the Mirage input catalogs. 
        
        # In[82]:
        
        
        contents = {}
        contents['1'] = {'times': times,
                         'fluxes': flux}
        
        
        # In[83]:
        
        
        # Save
        save_tso(contents, lightcurve_file, time_unit='second')
        
        
        # <a id="imaging_tso_catalog"></a>
        # ### Create Imaging TSO catalog
        
        # In this case, rather than the `tso_grism_catalog` in the yaml file, the user must supply a `tso_imaging_catalog`. This catalog will contain the list of sources whose flux will be varying with time. As in the grism TSO case, objects listed in the other types of catalogs will be added to the simulation in order to create a more realistic scene.
        # 
        # See the [notebook on catalog creation](https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb) in the examples directory of the Mirage repository for examples of how to create source catalogs. Catalogs may also be created manually.
        # 
        # Mirage does not yet support the use of input spectra in order to determine filter-based magnitudes for imaging simulations. Until that ability is added, users must estimate the source magnitude in each filter.
        
        # In[84]:
        
        
        imaging_tso_catalog = os.path.join(output_dir, 'tso_imaging_source.cat')
        
        
        # In[85]:
        
        
        tsimg_cat = ImagingTSOCatalog(ra=[object_ra], dec=[object_dec], lightcurve_file=[lightcurve_file])
        
        
        # In[86]:
        
        
        object_f182m_mag = 7.25
        object_f210m_mag = 7.23
        object_f470n_mag = 7.30
        
        
        # In[87]:
        
        
        # Add source magnitudes
        tsimg_cat.add_magnitude_column([object_f182m_mag], magnitude_system='vegamag',
                                       instrument='nircam', filter_name='f182m')
        tsimg_cat.add_magnitude_column([object_f210m_mag], magnitude_system='vegamag',
                                       instrument='nircam', filter_name='f210m')
        tsimg_cat.add_magnitude_column([object_f470n_mag], magnitude_system='vegamag',
                                       instrument='nircam', filter_name='f470n')
        
        
        # In[88]:
        
        
        tsimg_cat.save(imaging_tso_catalog)
        self.imaging_tso_catalog = imaging_tso_catalog

        
        # In[89]:
        
        
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
        
        # In[90]:
        
        
        # filter_list = ['F444W', 'F322W2','F182M','F210M','F470N']
        # box_width=140 ## arcsec
        # #                                                        besancon_catalog_file=besancon_result,
        
        # cat, mag_column_names = create_catalog.get_all_catalogs(object_ra, object_dec, box_width,
        #                                                         instrument='NIRCAM', filters=filter_list,
        #                                                         ra_column_name='RAJ2000', dec_column_name='DECJ2000',
        #                                                         starting_index=1)
        
        
        # In[91]:
        
    def prep_backgrounds(self):
        
        bbkgd_sources_ra = self.sys_params['background']['bbkgd_sources_ra']
        bkgd_sources_dec = self.sys_params['background']['bkgd_sources_dec']
        bkgd_sources_f444w_mag = self.sys_params['background']['bkgd_sources_f444w_mag']
        bkgd_sources_f322w2_mag = self.sys_params['background']['bkgd_sources_f322w2_mag']
        bkgd_sources_f182m_mag = self.sys_params['background']['bkgd_sources_f182m_mag']
        bkgd_sources_f210m_mag = self.sys_params['background']['bkgd_sources_f210m_mag']
        bkgd_sources_f470n_mag = self.sys_params['background']['bkgd_sources_f470n_mag']
        
        
        # In[92]:
        
        
        bkgd_cat_file = os.path.join(output_dir, 'ptsrcs.cat')
        
        
        # In[93]:
        
        
        bkgd_cat = PointSourceCatalog(ra=bkgd_sources_ra, dec=bkgd_sources_dec)
        
        
        # In[94]:
        
        
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
        
        
        # In[95]:
        
        
        bkgd_cat.save(bkgd_cat_file)
        
        
        # In[96]:
        
        
        bkgd_cat.table
        
        self.bkg_cat_file = bkg_cat_file

    def create_yaml(self):
        # <a id="yaml_files"></a>
        # ### Create input yaml files for Mirage
        
        # Mirage requires one input yaml file for each detector and exposure. Using the files exported from ATP combined with the various files produced above, Mirage can automatically generate all of the yaml files necessary to produce all of the data in the APT file. See the Mirage documentation on [creating yaml files from APT inputs](https://mirage-data-simulator.readthedocs.io/en/latest/yaml_generator.html) for details on the input parameters used below.
        
        # Populate the input catalog dictionary with the names of the source catalogs created above.
        
        # In[97]:
        
        
        catalogs = {'WASP-43': {'point_source': self.bkgd_cat_file,
                                'tso_imaging_catalog': self.imaging_tso_catalog,
                                'tso_grism_catalog': self.grism_tso_catalog,
                                }
                   }
        
        
        # Set the desired background level for the observations. Most common will to specify 'low', 'medium', or 'high'. The definitions of these three levels match those used in the ETC. Note that for the grism observations, if the `dateobs_for_background` is set to True, the background will be created based on the given observation date rather than the level given below.
        
        # In[98]:
        
        
        background = self.sys_params['background']['level']
        
        
        # Set the desired telescope roll angle and date for the observation. 
        
        # In[99]:
        
        
        pav3 = self.sys_params['obs']['pav3']
        dates = self.sys_params['obs']['dates']
        
        
        # In[100]:
        
        
        yam = yaml_generator.SimInput(self.xml_file, self.pointing_file, catalogs=catalogs, verbose=True,
                                      output_dir=self.output_yaml_dir, simdata_output_dir=self.output_data_dir,
                                      background=background, roll_angle=pav3,
                                      dates=dates, datatype='linear, raw', dateobs_for_background=True,
                                      reffile_defaults='crds')
        
        yam.use_linearized_darks = True
        yam.create_inputs()
        
        
        # List the yaml files created by the generator
        
        # In[101]:
        
        
        yam.yaml_files
        
        
        # # NOTE: I manually edited the date in the YAML file to
        # ```
        # date_obs: '2022-04-25'  # Date of observation
        # time_obs: '02:50:06.00'  # Time of observation
        # ```
        
        # In[106]:
        
        
        output_yaml_dir
        
        
        # In[111]:
    def create(self):
        
        orig_A5yaml_path = os.path.join(self.output_yaml_dir,'jw00042001001_01101_00001_nrca5.yaml')
        paramDict = yaml.safe_load(open(orig_A5yaml_path))
        
        
        # In[112]:
        
        
        paramDict['Output']['date_obs'] = self.sys_params['obs']['date_obs']
        paramDict['Output']['time_obs'] = self.sys_params['obs']['time_obs']
        
        
        # In[116]:
        
        
        #tmpSaveDir = os.path.join(output_yaml_dir,'test_yaml.yaml')
        yaml.safe_dump(paramDict,open(orig_A5yaml_path,'w'),default_flow_style=False)
        
        
        # Create a table showing details of what exposure each yaml file describes
        
        # In[117]:
        
        
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
        
        # In[118]:
        
        
        yaml.safe_dump(info_tab,open(os.path.join(self.output_dat_dir,'info_tab.yaml'),'w'),
                       default_flow_style=False)
        
        
        # <a id="create_simulated_data"></a>
        # ## Create Simulated Data
        
        # Simulate data from some of the files above.
        
        # <a id="grism_data"></a>
        # ### Grism TSO Data
        
        # First let's simulate a grism time series exposure using the F322W2 filter. 
        
        # In[119]:
        
        
        gr_tso_yaml_file = os.path.join(self.output_yaml_dir, 'jw00042001001_01101_00002_nrca5.yaml')
        
        
        gr_f322w2 = GrismTSO(gr_tso_yaml_file, SED_file=sed_file, SED_normalizing_catalog_column=None,
                            final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                            extrapolate_SED=True, override_dark=None, disp_seed_filename=None,
                            orders=["+1", "+2"])
        
        
        # In[ ]:
        
        
        # NOTE: This cell will take a while (> ~30 min) to run 
        gr_f322w2.create()
        
    def do_all():
        self.prep_sources()
        self.prep_background()
        self.create_yaml()
        self.create()
# <a id="accompanying_imaging_data"></a>
# ### Accompanying Imaging TSO Data

# Let's look at the imaging time series data collected on one of the shortwave detectors while the grism data above was being collected in the LW channel.

# In[64]:

#
#img_tso_sw_yaml = os.path.join(output_yaml_dir, 'jw88888001001_01101_00001_nrca1.yaml')
#
#
## In[ ]:
#
#
#img_tso = ImgSim()
#img_tso.paramfile = img_tso_sw_yaml
#img_tso.create()
#
#
## Examine the raw output file. This is the simulated data before any calibration. Detector A1 covers the same area as the [left half of detector A5](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview) (which was used to collect the Grism TSO above). This is why the source appears towards the right edge of A1. The +8 weak lens was also used in this observation, resulting in the expanded PSF.
#
## In[ ]:
#
#
## img_seed_file = os.path.join(output_data_dir,
##                              'jw88888001001_01101_00001_nrca1_uncal_F182M_seg001_part001_seed_image.fits')
#
#img_seed_file = os.path.join(output_data_dir,
#                             'jw88888001001_01101_00001_nrca1_uncal_F182M_WLP8_seed_image.fits')
#
#
## In[ ]:
#
#
#img_seed = fits.getdata(img_seed_file)
#
#
## In[ ]:
#
#
#show(img_seed[0, 4, :, :], 'Noiseless Seed Image', min=0, max=3000)
#
#
## Zoom in near PSF to compare from my estimate for WLP8 position
#
## In[ ]:
#
#
#plt.imshow(img_seed[0, 4, :, :],vmin=0,vmax=500)
#plt.xlim(1800,2048)
#plt.axvline(2048-64,color='red')
#
#
## Again, let's look at signal levels in all of the integrations, using the noiseless seed image. Since there is no noise and no background, we'll do a simple sum of the signal in a small box. Here we see that the depth of the transit matches that from the lightcurve we created with the Batman parameters above.
#
## In[ ]:
#
#
#signals = []
#cx, cy = 1857, 151
#for i in range(img_seed.shape[0]):
#    signals.append(np.sum(img_seed[i, 4, cy-1: cy+2, cx-1: cx+2]))
#
#
## In[ ]:
#
#
#f, a = plt.subplots()
#a.plot(np.arange(len(signals)), signals, color='red')
#a.set_xlabel('Integration Number')
#a.set_ylabel('Aperture Signal')
#
#
## Now look at the raw data
#
## In[ ]:
#
#
#uncal_file = os.path.join(output_data_dir, 'jw88888001001_01101_00001_nrca1_uncal.fits')
#
#
## In[ ]:
#
#
## head = fits.getheader(uncal_file)
## head
#
#
## In[ ]:
#
#
#data = fits.getdata(uncal_file)
#
#
## In[ ]:
#
#
#data.shape
#
#
## By displaying the difference between the last and first groups of one integration, we can remove bias signal that would otherwise make the object difficult to see.
#
## In[ ]:
#
#
#cds = 1.* data[0, 4, :, :] - 1. * data[0, 0, :, :]
#
#
## The image appears to have 4 512-column wide blocks because each block has been read out by a different amplifier. These differences will be minimized once the data have been run through the calibration pipeline. Note that this observation was taken with one of NIRCam's weak lenses in place, which is why the star's PSF has become so large.
#
## In[ ]:
#
#
#show(cds, 'Difference Image', min=0, max=1000)
#
#
## In[ ]:
#
#
#
#
