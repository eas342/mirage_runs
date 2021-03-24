# Mirage Runs
This has some TSO mirage simulations. You must install mirage first (https://mirage-data-simulator.readthedocs.io/en/latest/install.html).


You can grab this code like so
``` bash
git clone https://github.com/eas342/mirage_runs.git my_ERS_simulation
cd my_ERS_simulation
```

You can run a simulation for the ERS data challenge like so:

``` bash
nohup python NIRCam_ERS_full_02_new_trans.py > mirage_output.txt &
```
The nohup command allows it to run even if you are logged out. The output is also saved to `mirage_output.txt`.

