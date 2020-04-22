# HRR-BAM Data assimilation  

## Prerequisite
- Python3  
- [pyletkf](https://github.com/windsor718/pyletkf)  
- numpy  
- pandas  
- xarray (optional)  

## Quick start  
### Before you run a simulation  
1. Create configuration file by copying/editing config/config*.ini  
2. Compile Cython code for `dassim.py`  
```bash
python setup.py build_ext --inplace  
```
### Basics  
Run a new simulation:  
```bash
python dassim.py -c [configuration file] -s [start date] -e [end date]  -v
```
Run a new simulation without progressbar (recommended if you batch the script):  
```bash
python dassim.py -c [configuration file] -s [start date] -e [end date]  
```
Run a simulation from an existing restart file (restart.bin/restartAssim.bin in your output directory):  
```bash
python dassim.py -c [configuration file] -s [start date] -e [end date] --restart
```
Show help:
```bash
python dassim.py -h
```

### Change Log  
2020/03/02  
- data assimilation IO was updated. now it outputs reaches only in output_calibration.txt  
- add functions to read output_calibration.txt  
- add analysis/ensmean.py  
  
2020/02/27
- restarting IO scheme was updated in dassim.py. Native Numpy memory mapping was used for an efficiency. Now 1,000 times faster than before for this part.
- Overall computational efficiency was improved by more than 10 times. (e.g., for 30K reaches, total calc. time for the 10-year assimilation case is less than three hours with 20 cores)  
  
2020/02/26  
- Code was cleaned for an easier use.  
- Binary restart files IO was supported in HRR code. Now pyHRR is 3.5 ~ 4.0 times faster than before.  
