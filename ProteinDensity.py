# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 19:03:57 2022

@author: Dani
Using the same dictonary for creating the LAMMPS input
*.data file where the amino acids are in alphebtical order
Takes a FASTA file input and prints out the protein's molecular
weight, average residue weight, and space in vacuum. 
"""

import os #for adjusting system path parameters
import json #for parameter files
import sys #for adjusting command line arguments
import numpy as np
import Bio.SeqIO as seq

os.chdir(r"FILE LOCATION")
infile="*.FASTA"
    


#%% meantime initialization
AAs={}

# this cell only exists for now to establish default AAs parameters as those in old script
names=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
shortnames=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
charges=[0.,1.,0.,-1.,0.,0.,-1.,0.,.5,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.]
masses=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.17,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
sizes=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
lambdas=[.73,0.,.432,.378,.595,.514,.459,.649,.514,.973,.973,.514,.838,1.0,1.0,.595,.676,.946,.865,.892]
types=list(range(1,len(names)+1))
for i in range(len(names)):
    AAs[names[i]]={'charge':charges[i],'mass':masses[i],'type':types[i], 'size':sizes[i], 'lam':lambdas[i], 'shortname': shortnames[i]}
keys=list(AAs.keys())

#%%

sequence=list(seq.parse(infile,'fasta'))
reverse_lookup={}
for i in AAs.keys():
    reverse_lookup['{0}'.format(AAs[i]['shortname'])]=i
    
atomnum=0   
size=[]
mass=[]

for i in sequence[0].seq:
    ref=reverse_lookup[i]
    atomtype=AAs[ref]['type']
    size=np.append(size, AAs[ref]['size'])
    mass=np.append(mass, AAs[ref]['mass'])


MW = 0 
    
for i in range(len(mass)):
    MW = MW + mass[i]
    
Vol = 0

for i in range(len(size)):
    Vol = Vol + (4/3)*3.14*(size[i]/2)**3


MW_correct = MW # * (1/(6.022*10**20))
ProteinDensity = MW/Vol # * (1/(6.022*10**20)) # * (10**24)

averageaminoacidweight = MW / len(mass)

print("Total weight: ", MW_correct, 'Da and average amino acid weight' , averageaminoacidweight, 'Da')

#%%
# I need to determine "protein density"
# This is the density of the protein by it's self
# Size = diameter in angstroms
# Mass = molecular weight of each aminio acid in g/mol

Volume = 0
for i in range(len(size)):
    Volume = Volume + (4/3)*3.14*(size[i]/2)**3
    

Weight = 0   
for i in range(len(mass)):
    Weight = Weight + mass[i]

ProteinDensity = Weight/Volume
ProteinDensity_mgmL = ProteinDensity * (1/6.022E23) * 1000 * 1E24
print(ProteinDensity, "g/mol*A^3", ProteinDensity_mgmL, "mg/mL")

print(Volume)