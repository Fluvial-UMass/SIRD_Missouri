# pyHRR
A simple python wrapper to handle the HRR model. DAssim.py is a sample code to conduct data assimilation coupled with the pyHRR.   

hillslope_test: output_calibration is only contains validation reaches. -> pytHRR.py
- Note that this is still a test implementation. Data assimilation still stick with all reaches, so it dumps results for all reaches.  
- This will collapse the format of your output (some raws has only len(validation reaches) columns, whereas some raws has len(all reaches) columns).  
- Please fix this via analysis/ensProcess.py. This will remove non-validation reaches.  

