# SNHostCheck
Quick script to cross match young SN candidates with PS1 stacked photometry for hosts. 
Assumes something with PS1 photometry within a radius od 0.1"<r<2.0" is a host.

This script requires python 3 and astropy, alongside other classical libraries like numpy and pandas.

You must have an output .csv from the YSE-PZ SQL Explorer: https://ziggy.ucolick.org/yse/explorer/202/

To run:

    python SNHostCheck.py New_young_objects.csv
