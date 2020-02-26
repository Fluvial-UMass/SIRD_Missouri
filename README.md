# HRR-BAM Data assimilation [Internal]  

## Prerequisite  
- [pyletkf](https://github.com/windsor718/pyletkf)  

## Quick start  
### Before you run a simulation  
1. Create configuration file by copying/editing config/config*.ini  
2. Compile Cython code for `dassim.py`  
```python
python setup.py build_ext --inplace  
```
### Basics  
Run a new simulation:  
```python
python dassim.py -c [configuration file] -s [start date] -e [end date]  
```
Run a simulation from existing restart file (restart.txt/restartAssim.txt in your output directory):  
```python
python dassim.py -c [configuration file] -s [start date] -e [end date] --restart
```
Show help:
```python
python dassim.py -h
```