#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import numpy as np
from pyopenms import *
#_________Modify the file handler location 
fh = open("y1.fasta")
#get the sequence from file
bsa = "".join([l.strip() for l in fh.readlines()[1:]])
# getting contamients from file
tr_fh = open("tr.fasta")
trypsine = "".join([l.strip() for l in tr_fh.readlines()[1:]])
#joining values into strings
bsa = AASequence.fromString(bsa+trypsine)
#digest  the yeast sequence 
digsted_peptides = []
dig = ProteaseDigestion()
dig.digest(bsa, digsted_peptides,3,20)
print("peptides number : " ,len(digsted_peptides))

#getting experimental spec
exp = MSExperiment()
MzMLFile().load("xp45.mzML", exp)#exp file
exp_spec = exp.getSpectra()
first_exp_spec = exp_spec[0]

#preparing the theortical mass spec
spec = MSSpectrum()
tsg = TheoreticalSpectrumGenerator()
p = Param()
p.setValue("add_b_ions", "false")
p.setValue("add_y_ions", "true")
p.setValue("add_a_ions", "false")
p.setValue("add_losses", "false")
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
for peptide in digsted_peptides[0:]:  # <_____________________________TO control number of peptides 
    tsg.getSpectrum(spec, peptide, 1, 2) 
#for ion, peak in zip(spec.getStringDataArrays()[0], spec):
    #1+1
    #print(ion.decode(), "is generated at m/z", peak.getMZ())    
plt.xlabel("m/z")
plt.ylabel("intensity")
plt.bar(spec.get_peaks()[0], spec.get_peaks()[1])

#drawing mirror plot function 
def mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title):
    obs_int = [element / max(obs_int) for element in obs_int] # relative intenstiy
    theo_int = [element * -1 for element in theo_int] # invert the intensity for the mirror plot
    plt.figure(figsize=(15,10))
    plt.bar(obs_mz, obs_int, width = 3.0)
    plt.bar(theo_mz, theo_int, width = 3.0)
    plt.title(title)
    plt.ylabel('intensity')
    plt.xlabel('m/z')
#initialzing the experimental  mz & intinsite
exp_mz, exp_int = first_exp_spec.get_peaks()
print("___________________________________\n")  
min_range=min(exp_mz) 
max_range=max(exp_mz)
# We filter the peaks of the theoretical spectrum to fit the range (to reduce image complexity) by 5%
max_range= max_range -(max_range*0.05)
min_range= min_range -(min_range*0.05)
theo_mz, theo_int = [], []
for mz, intensity in zip(*spec.get_peaks()):
    if mz >= min_range and mz <= max_range:
        theo_mz.append(mz)
        theo_int.append(intensity)

title = 'experimental vs theoretical spectrum'
mirror_plot(exp_mz, exp_int, theo_mz, theo_int, title)
                                    # the matching alignment
alignment = []
spa = SpectrumAlignment()
ppp = spa.getParameters()
ppp.setValue("tolerance", 0.5)#using default 0.5 dalton tolerance
spa.setParameters(ppp)
# align both spectra
sortedspec=spec.sortByPosition()
spa.getSpectrumAlignment(alignment, spec, first_exp_spec)
# Print matching ions and mz from theoretical spectrum
print("Number of matched peaks: " + str(len(alignment)))
print("___________________________________\n")  
print("ion\ttheo. m/z\tobserved m/z")
for theo_idx, obs_idx in alignment:
    ion_name = spec.getStringDataArrays()[0][theo_idx].decode()
    ion_charge = spec.getIntegerDataArrays()[0][theo_idx]
    print(ion_name + "\t" + str(ion_charge) + "\t"
      + str(spec[theo_idx].getMZ())
      + "\t" + str(first_exp_spec[obs_idx].getMZ()))
    

    
theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
for theo_idx, obs_idx in alignment:
    theo_mz.append(spec[theo_idx].getMZ())
    theo_int.append(spec[theo_idx].getIntensity())
    obs_mz.append(first_exp_spec[obs_idx].getMZ())
    obs_int.append(first_exp_spec[obs_idx].getIntensity())
title = 'Observed vs theoretical spectrum (aligned)'
mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)

plt.show()
# In[ ]:




