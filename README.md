# UDSSNe
Analysis of the UDS dataset to search for supernovae at high redshifts

.gitignore           excludes large, produced datafiles from being uploaded
README.md            this file

detection.py         Chi squared test applied to UDS dataset in both year an monthly data. Can be modified for any dataset to run chi squared on.
import.py            imports both monthly and yearly data from UDS
sample generator.py  generates samples to test algorithms on based off the UDS data
uds_data_handling.py alternate method of import.py. May not be consistent with file format.

data/                current working directory for large data products
Test Light Curve Library/   contains a bunch of sample lightcurves that may be used for the testing of our algorithm (maybe not likely)
