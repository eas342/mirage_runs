# Mirage Runs
This repository contains a wrapper to run JWST NIRCam simulations using the mirage tool. You must install mirage first (https://mirage-data-simulator.readthedocs.io/en/latest/install.html).

## Installation

You can install the prerequisites and exact versions of tools as follows. This assumes you already have a conda environment set up for mirage (in this case it is called `mirageRuns`), following the mirage instructions.

``` bash
git clone https://github.com/eas342/mirage_runs.git
cd mirage_runs ## or wherever you have placed this repository
conda activate mirageRuns
pip install -r requirements.txt

```
### Testing the installation on a minimal example.
After your installation is complete, you can test

``` bash
conda activate mirageRuns # if relevant environment is not already active
cd mirage_runs ## or wherever you have placed this git repository
nohup python NIRCam_grism_tso_wrapper.py > mirage_test_output.txt &
```

## Usage

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
4. Export a `.pointing` file using File -> Export. Make sure it is the same name as the `.xml` file. For example, `my_program.xml` and `my_program.pointing`
5. Create a parameter file as shown below.
6. TBD: Run the Mirage Wrapper

### Parameter File

### Specific code for the ERS data hackathon:

You can run a simulation for the ERS data challenge like so:

``` bash
nohup python NIRCam_ERS_full_02_new_trans.py > mirage_output.txt &
```
The nohup command allows it to run even if you are logged out. The output is also saved to `mirage_output.txt`.

