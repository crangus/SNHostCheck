#Candidate SN Host Info

# You need to run the SQL Explorer Script "" and download the csv file this produces.
#

from __future__ import print_function
import sys
sys.path.insert(0, './ps1/')
import argparse
import numpy as np
import scipy as sp
import pandas as pd
from ps1 import PS1_query
import astropy
import astropy.units as u

parser = argparse.ArgumentParser(description='Young Candidates:')
parser.add_argument('file', type=str, help='input csv file',default=None)
condits=parser.parse_args()
file = condits.file

candidates=pd.read_csv(file,sep=',',header=0)

print("Ingesting %i recent candidates from YSE-PZ"%(len(candidates.name)))
no_host_flag=[]
host_flag=[]

print('Checking host environments')
for row in np.arange(0,len(candidates),1):
    nearby=PS1_query.search(candidates.name.iloc[row],candidates.ra.iloc[row],candidates.dec.iloc[row],1).to_pandas().sort_values('Seperation').reset_index()
    if len(nearby)==0:
        no_host_flag.append(candidates.name.iloc[row])
    if len(nearby)>0:
        if nearby.iloc[0]['Seperation']>0.1:
            if nearby.iloc[0]['Seperation']<=2:
                    host_flag.append(candidates.name.iloc[row])

print('Total objects meeting all selection criteria (%i/%i):'%(len(host_flag),(len(candidates.name))))

for i in host_flag:
	print(i)
