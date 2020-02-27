# HRR-BAM Data assimilation [Internal]  

## Prerequisite
- Python3  
- [pyletkf](https://github.com/windsor718/pyletkf)  
- numpy  
- pandas  

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
Run a simulation from an existing restart file (restart.txt/restartAssim.txt in your output directory):  
```bash
python dassim.py -c [configuration file] -s [start date] -e [end date] --restart
```
Show help:
```bash
python dassim.py -h
```

### Change Log  
2020/02/27
- restarting IO scheme was updated in dassim.py. Native Numpy memory mapping was used for efficiency. Now 1,000 times faster than before.
- Overall computational efficiency was improved by 40 times. (e.g., for 30K reaches, total calc. time for the assimilation case is less than a hour with 20 cores)  
2020/02/26  
- Code was cleaned for an easier use.  
- Binary restart files IO was supported in HRR code. Now pyHRR is 3.5 ~ 4.0 times faster than before.  
