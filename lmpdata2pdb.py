# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 13:39:58 2024

@author: Dani
The following script takes LAMMPS *.data files and turns them into *.pdb files
This works for coarse-grained simulations where each amino acid is represented by a single bead
The center of mass for the 'bead' is replaced with an alpha carbon

Only *.data, the number of chains, and residue identifier numbers are required
The output file is simply a *.pdb
"""

import os #for adjusting system path parameters
import json #for parameter files
import sys #for adjusting command line arguments
import numpy as np
import Bio.PDB as pdb
import pandas as pd
import string

#%%
# Here you place the information required for input for the script 
os.chdir(r"FILE LOCATION") # in final iteration this file will be run from the command line, so it will be run in the working directory anyway
data_file = '*.data'
number_of_chains = 100
output_file = '*.pdb'

#%%
AAs={}
# dictionary to convert LAMMPS ID number into residue type for *.pdb
# Please ensure that the names, masses, and sizes match.
# the pre-made dictonary has the residues in alphebtical order. 
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

counts = (list(range(1, len(names))))
AA_name = dict(zip(counts,names))

# Function to generate alphabetical labels for chain IDs for *. above 62 chains
def generate_labels(n):
     chainID = []
     for i in range(n):
         label = ''
         while i >= 0:
             label = string.ascii_uppercase[i % 26] + label
             i //= 26
             i -= 1
         chainID.append(label)
     return chainID


# Number of elements = number of chains
num_elements = number_of_chains
# Generate labels for *.pdb Chain IDs
chain_identifiers = list(string.ascii_uppercase) + list(string.digits) + list(string.ascii_lowercase) + generate_labels(num_elements)[26:]
ChainID = chain_identifiers[:num_elements+1]
counts = (list(range(1, num_elements+1)))
Chain_dic = dict(zip(counts, ChainID))

#loads the *.data into a pandas dataframe for easy manipulation 
df = pd.read_csv(data_file, sep='\ ', skiprows=35, names=['atom#','chain#','AA#','charge','x','y','z','a','b','c'])
df = df.dropna() #drops nan values given by bond information within data file
df['atom#'] = df['atom#'].astype(int) #makes the atom# integers instead of objects for sorting
df_sorted = df.sort_values(by=['atom#'], ascending=[True]) #puts atoms# in correct order for *.pdb

df_sorted['AA#']=df_sorted['AA#'].replace(AA_name) #gives residues correct 3 digit naming for *.pdb from 1 number coarse grain ID
df_sorted['chain#']=df_sorted['chain#'].replace(Chain_dic) #replaces chain ID from *.data with correct naming for *.pdb
df_sorted = df_sorted.reset_index(drop=True) #reorders index to ensure atom# does not reset
df_sorted['atom_type'] = 'CA' #adds alpha carbon atom type for all beads


#%%
# Below writes the output file ensuring correct spacing for *.pdb
file = open(output_file,'w')
k = 1
for i in range(len(df_sorted)-1):
    print('ATOM  ', f"{df_sorted.at[i,'atom#']:>5} {'CA':<4}{df_sorted.at[i,'AA#']:>3}{df_sorted.at[i,'chain#']:>2} {k:>4} {df_sorted.at[i,'x']:>8.3f}{df_sorted.at[i,'y']:>8.3f}{df_sorted.at[i,'z']:>8.3f}", file=file)
    k = k+1
    if df_sorted.at[i,'chain#'] != df_sorted.at[i+1,'chain#']: #this is why the last chain is not picked up... 
       print("TER", file=file)
       k = 1
k = k+1
print('ATOM  ', f"{df_sorted.at[len(df_sorted)-1,'atom#']:>5} {'CA':<4}{df_sorted.at[len(df_sorted)-1,'AA#']:>3}{df_sorted.at[len(df_sorted)-1,'chain#']:>2} {k:>4} {df_sorted.at[len(df_sorted)-1,'x']:>8.3f}{df_sorted.at[len(df_sorted)-1,'y']:>8.3f}{df_sorted.at[len(df_sorted)-1,'z']:>8.3f}", file=file)
print("TER", file=file)
file.close()
