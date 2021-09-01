# Mirage Runs
This repository contains a wrapper to run JWST NIRCam simulations using the mirage tool. You must install mirage first (https://mirage-data-simulator.readthedocs.io/en/latest/install.html). The intent of this repository is to avoid the version control problem of using the TSO notebook and tracking down/debugging what was different about each mirage run. It also puts all of the input parameters in one place to separate it from the code of creating the mirage simulation.

## Installation

You can install the prerequisites and exact versions of tools as follows. This assumes you already have a conda environment set up for mirage (in this case it is called `mirageRuns`), following the mirage instructions.

``` bash
git clone https://github.com/eas342/mirage_runs.git
cd mirage_runs ## or wherever you have placed this repository
conda activate mirageRuns
pip install -r requirements.txt

```
### Testing the installation on a minimal example.
After your installation is complete, you can test it. This should take a few minutes to finish a simple simulation of a few short integrations.

``` bash
conda activate mirageRuns # if relevant environment is not already active
cd mirage_runs ## or wherever you have placed this git repository
nohup python NIRCam_grism_tso_wrapper.py > mirage_test_output.txt &
```

Check `mirage_test_output.txt` for any errors and if not, whether it says `Grism TSO simulator complete` at the bottom.
Also, check to see if `mirage_output_001/sim_data` contains grism file named `jw88888001001_01101_00002_nrca5_uncal.fits`.

## Step-by-Step instructions

There are a handful of steps needed to run a simulation. Be prepared to gather

* detailed exposure information on your observation, such as from the <a href="https://jwst.etc.stsci.edu">Exposure Time Calculator (ETC)</a>.
* <a href="https://www.stsci.edu/scientific-community/software/astronomers-proposal-tool-apt">Astronomer's Proposal Tool (APT)</a>
* Planet parameters, including 
	+ stellar parameters
		- Teff
		- metallitiy
	   - logg
	+ the system's brighness
	   - K magnitude (Vega system)
	   - F444W magnitude (likely the same as the K mag within errors)
	   - F322W2 magnitude (likely the same as the K mag to within errors)
	+ the system's RA & Dec coordinates in degrees
	+ planet parameters
	   - transit center time
	   - period
	   - planet radius/stellar radius Rp/R*
	   - semi-major axis/stellar radius, a
	   - inclination in degrees
	   - eccentricity
	   - longitude of periastron (if eccentricity is not 0)
	   - limb darkening law
	   - A theoretical planet spectrum (or one can assume a constant radius)
   + background star information
   	   - coordinates Ra & Dec in degrees
   	   - Brightness in NIRCam bands
   + observation dates


Then, follow these stps

1. Calculate the exposure parameters using the <a href="https://jwst.etc.stsci.edu">ETC</a> and fill out the <a href="https://www.stsci.edu/scientific-community/software/astronomers-proposal-tool-apt">APT</a> file.
2. Start APT and create your exposures for the NIRCam grism time series template.
3. Export a `.xml` file using File -> Export
4. Export a `.pointing` file using File -> Export. **Make sure it is the same name as the `.xml` file**. For example, `my_program.xml` and `my_program.pointing`
5. Create a parameter file as described below.
6. Run the Mirage Wrapper as a script. For example:

``` bash
conda activate mirageRuns # if relevant environment is not already active
cd mirage_runs ## or wherever you have placed this git repository
nohup python NIRCam_grism_tso_wrapper.py input/wasp43.yaml > output/mirage_out_wasp43.txt &
```

## Parameter File

You should create a directory, ie `mkdir mirage_input_003_wasp43` and put a `.yaml` file in it, following the example: <a href='mirage_input_001/source_params.yaml'> `mirage_input_001/source_params.yaml`</a>.

### descrip
A string to describe the mirage Run. For example `"WASP-43 b Transit"`

### stellar

`teff`: The stellar effective temperature (K). For example 5600.

`metallicity`: The stellar metallicity ([Fe/H]). For example 0.0 for solar.

`logg`: The stellar surface gravity log_10(gravity/ (g/cm^2)). For example, 4.5.

### system
`kmag`: The 2MASS Vega Magnitude of the system in the K_S band (the default for looking something up in Simbad.

`ra`: The right ascension coordinates in **degrees**. For example, `66.3709122`

`dec`: The declination coordinates in **degrees**. For example `-30.6000045`.

`f444w_mag`: The Vega magnitude in the JWST F444W band. In most cases, you can set this equal to the `kmag` because it is the Raleigh Jeans limit for stars and is probably the same within errors.

`f322ww_mag`: The Vega magnitude in the JWST F444W band. In most cases, you can set this equal to the `kmag` because it is the Raleigh Jeans limit for stars and is probably the same within errors.

### planet
`t0`: The epoch at transit center **relative to the start of the exposure**, in seconds.

`per`: The orbital period, in seconds. For example 3600. * 24. * 3.0 = 259200.

`rp`: The ratio of the planet radius to stellar radius. For example, 0.1.

`a`: The ratio of the semi-major axis to stellar radii, a/R*. For example, 9.5.

`inc`: The inclination of the planet orbit, in degres. For example, 89.1.

`ecc`: The eccentricy of the planet orbit. For example, 0.02

`w`: The longitude of pariastron for the planet orbit. If eccentricty is 0.0, set this so 90.0

`limb_model`: The limb darkening model. For example, `"quadratic"` or `"nonlinear"`

`ld_u`: The limb darkening values. For example, `[0.1,0.05]` for `"quadratic"` or `[0.1,0.1,0.1,-0.05]` for `"nonlinear"`

`spec_file`: The path to the spectrum file that contains the planet spectrum.

### background
This describes both the zodiacal background as well as any other neighboring stars.

`level`: (string) The level of the Zodiacal background. For example, `"medium"`

`bkgd_sources_ra`: A list of right ascensions of background sources **in degrees**. For example, `[66.370916666,66.380916666]`.

`bkgd_sources_dec`: A list of delclinations of background sources **in degrees**. For example, `[-30.5994003, -30.5994009]`.

`bkgd_sources_f444w_mag`: Magnitude of background source in the F444W band. And so and so-forth for `f3222w2_mag`, `f182m_mag`, `f210m_mag`, `f470n_mag`.

### obs
Describes the observations.

`dur`: (float) Duration of the exposure used to calculate interpolating functions.

`pav3`: (float) The position angle of the V3 axis at observation in degrees. This is important when there are neighboring stars that could contaminate the spectrum of the target.

`date_obs`: (string) The exact date of the exposure start. For example, `"2022-11-01"`.

`time_obs`: (string) The exact time of the exposure start. For example, `"02:50:06.000"`. Note that the transit center will happen at this time plus the value above under `planet: t0`.

`obsnum`: (list of ints) The observation numbers from the APT file. If "null" in the yaml file (`None` in Python), it will use the first one available. For example, `[1]`.

### dataPaths
These say where to grab input data and where the outputs are.

`input_data_path`: (string). Directory for input data: XML, source parameters. This
            is also where Mirage will dump mirage yaml files

`output_data_dir`: (string).            Directory for output data where the simulations and intermediate files will go.

`xmlFile`: (string). Name of the XML file made by APT. Put that and `.pointings` in the input path

## Specific code for the ERS data hackathon:

There is a script for a specific (old) mirage simulation for an <a href="https://ers-transit.github.io/pre-launch-hackathon.html">Early Release Science (ERS) data hackathon</a>. This uses an older version of mirage and inputs. You can run a simulation for the ERS data challenge like so:

``` bash
nohup python NIRCam_ERS_full_02_new_trans.py > mirage_output.txt &
```
The nohup command allows it to run even if you are logged out. The output is also saved to `mirage_output.txt`. 

