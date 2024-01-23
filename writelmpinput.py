# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 17:03:09 2020

@author: glham with slight modifications from Dani

rewrote the original script to more generally handle pdb files without neading to remove headers/make sure only one chain is present, etc
Decided to use biopython in this version for 3 main reasons:
    1) consistency in opening PDB files is built right in, and extensible to other formats as generated structure objects are compatible with them
    2) makes parameters easier to call than manually parsing PDB files
    3) Streamline process of simply copying parameters from initial file for weird molecules/atoms
    4) not a real reason, just wanted to practice it
"""
#biopython
import os #for adjusting system path parameters
import json #for parameter files
import sys #for adjusting command line arguments
import numpy as np
import Bio.PDB as pdb
#%%
#import conversionfunctions as cf

os.chdir(r"FILE LOCATION") # in final iteration this file will be run from the command line, so it will be run in the working directory anyway

# initialize parameters
# initially, here is where parameters will be defined. in the command line version, this cell will
# instead pull from the command line arguments with tags like -i (input) -o (output) 
#initial script will be written adduming we want 1 bead model with beads located at the Ca positions for each residue
#future options: dif jsons for parameters sets, either 1 chain or all chains, include bonds, angles, and dihedrals
#also want an option to build a model of a linear chain (at first only 1 bead, as setting up the proper bond angles will be a nightmare)
name = "FileName"

infile=name+".FASTA"#3gsl.pdb' #full path relative to working directory
savemodel=False
#outmodel='3gslbead' #the output altered pdb file. not supported yet, will implement cell that builds the PDB as either a struct object or just writes a file
saveLAMMPin=True
outfile=name+'_lmpinput' #the output datainfile for LAMMPS
lmpout=name+'_parameters' #input script for LAMMPS with pair coeffs, etc (currently just pair coeffs to copy past into example)
informat='fasta' #'pdb'
num_chain=1
duplicates=100 #for 1 chain, number of copies to put in box. for multiple, sampples the dif. chains to make duplicates (not implemented). translates duplicates around in box
remove_nonAAs=True #recommended, since most likely ions, etc. can just be added back later. Also, treating non-AAs is currently not supported for things not in the AAs dictionary
edit_atom_properties=True #whether to alter atomic properties using Atoms dictionary (currently called AAs, should change for generality and avoid confusion "amino acids" vs "atom attributes")
write_bonds=True
write_angles=False
write_dihedrals=False
box_xyz_dimensions=[1000.,1000.,1000.] # half total width of box in each dimension
initial_position=-box_xyz_dimensions[0]/2 # for new chains, where to place the first atom in the box
make_ljlambda_coeffs=True
makexyz=True
placement_method='standard_size'
standard_size=3.6
shiftx=True

AAs={} #json/dictionary containing the parameters for "". Note this approach allows extension to unnatural amino acids simply by adding a new dictionary element
#Atoms={} #atomic properties when the input file properties are not to be preserved or don't exist
#Bonds={} #additional bonds beyond the defaults (defaults are linear chain for 1 bead, preserved atomic bonds for all atom, etc.)

conversion_type='1_bead_Ca' #future hope is multiple, including 1_bead_CoM, multiple_bead (such as Ca+sidechain), all atom, etc.

# meantime initialization
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


# Atoms_section write
# conversion to new model for atoms section
#here we'll just write the conversion script for conversion to one bead, but eventually the hope is to have a list of functions
# in a separate file (for cleanliness), and simply import them to use here in a single line before writing.

#This script version will be 1-bead CA based on the first chain in a pdb file
#after conversion, duplication will be another cell
Atoms_section=[]

if informat=='pdb':
    import Bio.PDB as pdb #for reading pdb files
    parser=pdb.PDBParser()
    structure=parser.get_structure("struc",infile)
    atomnum=0
    #assumes num_chain=1
    for residue in structure[0].child_list[0]:
        atomnum+=1
        if not pdb.is_aa(residue):
            if remove_nonAAs: #will change this to handle non-AAs later
                continue
        if conversion_type=='1_bead_Ca':
            loc=residue["CA"].coord
            chainnum=1 #when looping over multiple chains, want to grab this value
            atomtype=AAs[residue.resname]['type']
            atomcharge=AAs[residue.resname]['charge']
        elif conversion_type=='1_bead_CoM':
            continue #not implemented yet
        Atoms_section.append('\t {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \n'.format(atomnum,chainnum,atomtype,atomcharge,loc[0],loc[1],loc[2]))

if informat=='fasta': #only supports 1 bead model for the moment, and reads out first chain in sequence only
    import Bio.SeqIO as seq
    sequence=list(seq.parse(infile,informat))
    reverse_lookup={}
    for i in AAs.keys():
        reverse_lookup['{0}'.format(AAs[i]['shortname'])]=i
    atomnum=0
    chainnum=1
    size=[0.,0.] #store the atomsize for the current and previous atom
    loc=[initial_position,0.,0.] #position to place the atom at (this script adds them along x axis in a line)
    if duplicates>1:
        chain1p=[]
    if makexyz==True:
        xyz=[]
    for i in sequence[0].seq: #for now just use a single chain, the first in the list
        ref=reverse_lookup[i]   
        atomnum+=1
        atomtype=AAs[ref]['type']
        atomcharge=AAs[ref]['charge']
        size[0]=size[1]
        size[1]=AAs[ref]['size']
        if placement_method=='mean_size':
            loc[0]+=(size[0]+size[1])/2. #place the next atom a distance away from the previous equal to the average size of the two atoms
        elif placement_method=='standard_size':
            loc[0]+=standard_size
        if duplicates>1:
            chain1p.append([atomnum,chainnum,atomtype,atomcharge,loc[0],loc[1],loc[2]])
        if makexyz==True:
            xyz.append([loc[0],loc[1],loc[2]])
        Atoms_section.append('\t {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \n'.format(atomnum,chainnum,atomtype,atomcharge,loc[0],loc[1],loc[2]))


# Duplicate section     
#really badl programmed atm, just to get it done. ideally ill be incorporated above and with a better way than doing the above loop duplicate num times (basically, copy it instead after making Atoms_section)
chain1atomnum=atomnum
xlength=chain1p[-1][4]-chain1p[0][4]   
if duplicates>1.:
    for i in range(1,duplicates):
        chainnum+=1
        newy=np.random.uniform(-box_xyz_dimensions[1],box_xyz_dimensions[1])/2
        newz=np.random.uniform(-box_xyz_dimensions[2],box_xyz_dimensions[2])/2
        newxadd=np.random.uniform(0,box_xyz_dimensions[0]-xlength)
        for j in range(chain1atomnum):
            if shiftx:
                newx=chain1p[j][4]+newxadd
            else:
                newx=chain1p[j][4]
            atomnum+=1
            Atoms_section.append('\t {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \n'.format(atomnum,chainnum,chain1p[j][2],chain1p[j][3],newx,newy,newz))
            if makexyz==True:
                xyz.append([newx,newy,newz])
        

if makexyz==True:
    xyzfile=open(outfile+'.xyz',"w")
    xyzfile.write('{0} \n \n'.format(atomnum))
    for i in xyz:
        xyzfile.write('C \t {0} \t {1} \t {2} \n'.format(i[0],i[1],i[2]))
    xyzfile.close()


# header section write

    
#some of the following lines need to change as things get more complex, such as the number of bonds=/= atomnum-1 for more complex bond types
Header_section=[]
Header_section.append('\t {0} atoms \n'.format(chain1atomnum*duplicates))
if write_bonds:
    Header_section.append('\t {0} bonds \n'.format((chain1atomnum-1)*duplicates))
if write_angles:
    Header_section.append('\t {0} angles \n'.format((chain1atomnum-2)*duplicates))
if write_dihedrals:
    Header_section.append('\t {0} dihedrals \n'.format((chain1atomnum-3)*duplicates))

Header_section.append('\n')

Header_section.append('\t \t {0} atom types \n'.format(len(AAs.keys())))
if write_bonds:
    Header_section.append('\t \t {0} bond types \n'.format(1))
if write_angles:
    Header_section.append('\t \t {0} angle types \n'.format(1))
if write_dihedrals:
    Header_section.append('\t \t {0} dihedral types \n'.format(1))
    
Header_section.append('\n')

Header_section.append('\t -{0} {0} xlo xhi \n'.format(box_xyz_dimensions[0]))
Header_section.append('\t -{0} {0} ylo yhi \n'.format(box_xyz_dimensions[1]))
Header_section.append('\t -{0} {0} zlo zhi \n'.format(box_xyz_dimensions[2]))   

Header_section.append('\n')
    

# Masses section write
Masses_section=[]
    
for i in AAs.keys():
    Masses_section.append('\t {0} \t {1} \n'.format(AAs[i]['type'], AAs[i]['mass']))

# Bonds section write
track=0
if write_bonds:
    Bonds_section=[]
    for i in range(1,atomnum):
        if (i)%chain1atomnum==0:
            track+=1
            continue
        Bonds_section.append('\t {0} \t {1} \t {2} \t {3} \n'.format(i-track, 1, i, i+1))
    Bonds_section.append('\n')    
# Angles section write
track=0
if write_angles:
    Angles_section=[]
    for i in range(1,atomnum-1):
        if (i+1)%chain1atomnum==0 or (i)%chain1atomnum==0:
            track+=1
            continue
        Angles_section.append('\t {0} \t {1} \t {2} \t {3} \t {4} \n'.format(i-track, 1, i, i+1, i+2))   
    Angles_section.append('\n')    
# Dihedrals section write
track=0
if write_dihedrals:
    Dihedrals_section=[]
    for i in range(1,atomnum-2):
        if (i+1)%chain1atomnum==0 or (i+2)%chain1atomnum==0 or (i)%chain1atomnum==0:
            track+=1
            continue
        Dihedrals_section.append('\t {0} \t {1} \t {2} \t {3} \t {4} \t {5} \n'.format(i-track, 1, i, i+1, i+2, i+3))    
    Dihedrals_section.append('\n')    
    
# write in file

Lfile=open(outfile,"w")
Lfile.write('#HEADER \n \n')
for i in Header_section:
    Lfile.write(i)
Lfile.write('\nMasses \n \n')
for i in Masses_section:
    Lfile.write(i)
Lfile.write('\nAtoms \n \n')
for i in Atoms_section:
    Lfile.write(i)
if write_bonds:
    Lfile.write('\n \nBonds \n \n')
    for i in Bonds_section:
        Lfile.write(i)
if write_angles: 
    Lfile.write('\n \nAngles \n \n')   
    for i in Angles_section:
        Lfile.write(i)
if write_dihedrals:
    Lfile.write('\n \nDihedrals \n \n')
    for i in Dihedrals_section:
        Lfile.write(i)
Lfile.close()
# pair coefficients
#make pair coefficients file for the ljlambda FF to be copy/pasted into example file (*.lmp)
if make_ljlambda_coeffs==True:
    lmpfile=open(lmpout, "w")
    lmpfile.write('pair_coeff \t * \t * \t 0.000000 \t 0.000 \t 0.000000 \t 0.000 \t 0.000 \n')
    for i in range(len(keys)):
        for j in range(i,len(keys)):
            if AAs[keys[i]]['charge']!=0. and AAs[keys[j]]['charge']!=0:
                cutoff=35.000
            else:
                cutoff=0.000
            lmpfile.write('pair_coeff \t {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \n'.format(i+1,j+1,0.200000,(AAs[keys[i]]['size']+AAs[keys[j]]['size'])/2.,(AAs[keys[i]]['lam']+AAs[keys[j]]['lam'])/2.,(AAs[keys[i]]['size']+AAs[keys[j]]['size'])/2.*3.,cutoff))
    lmpfile.close()