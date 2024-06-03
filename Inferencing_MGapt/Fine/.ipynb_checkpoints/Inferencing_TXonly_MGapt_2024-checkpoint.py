#!/usr/bin/env python
# coding: utf-8

# # A chemical reaction network model of PURE
# Zoila Jurado [1,*], Ayush Pandey [1], Richard M. Murray [1,2]
# 
# August 8, 2023
# 
# [1] Division of Engineering and Applied Science, California Institute of Technology, Pasadena, CA \
# [2] Division of Biology and Biological Engineering, California Institute of Technology, Pasadena, CA \
# [*] Corresponding author: zjuradoq@caltech.edu

# In this notebook, we built a chemical reaction network transcription model for PURE protein synthesis. The transcription and translation models are separately described and generated. All species and reactions in the transcription and translation model are then combined together with additional reactions and species. The combined PURE model can be used to model the expression of a protein in PURE. We used this notebook to model the transcription and translation of a plasmid containing the construct: pT7-MGapt-UTR1-deGFP-tT7.
# 
# A full description of the model can be found in the BioRvix. 

# # Run inference

# In[4]:


"""To run this file on HPC, simply run:
`sbatch inference_run.sh`
The results are written to `mcmc_results.csv`
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import csv

from bioscrape.types import Model
from bioscrape.inference import py_inference
from bioscrape.simulator import py_simulate_model

#Get directory
import os
directory = os.getcwd()

# # Bayesian Inference 
from bioscrape.inference import py_inference

import matplotlib.pyplot as plt
import numpy as np
from bioscrape.types import Model
import time
import pandas as pd
import warnings
from scipy.integrate import odeint
import seaborn as sns
import pandas as pd


# In[6]:


# ## Getting experimental data
filename='CRN_PURE_MGaptTXonlyCourse'
file = '/Data_files/MGapt_mRNA_data_final.csv'
mGapt_Data =  pd.read_csv(directory+ file, delimiter = '\,', names = ['time','mRNA0', 'mRNA1','mRNA2'], engine='python',skiprows = 1)
first_row = mGapt_Data.iloc[0, 1:].mean()  # extract the first row, excluding the 'time' column
mGapt_Data.iloc[:, 1:] -= first_row  # subtract the values in the first row from all other rows, excluding the 'time' column
# Truncate at 4 hours
truncated_data = mGapt_Data[mGapt_Data['time'] <= 3*3600]
# Define timepoints
timepoints =  np.array(truncated_data['time'])

# Run model to test
m =Model(sbml_filename = filename+"_updated.xml")
initial_con={'ATP': 2206.0,'CTP': 1250.0, 'DNA': 0.005, 'T7RNAP':1,
             'GTP': 590.0,'UTP': 1250.0,}
m.set_species(initial_con)

# Simulate model
R0 = py_simulate_model(Model = m, timepoints = timepoints)
plt.figure(figsize=(6, 8))
plt.title("mRNA Concentration Over Time")
plt.plot(R0['time']/3600, R0['MGapt'], label = "mRNA", color = 'm') #mRNA

#Experimental data set 1
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA0'], label = "T7-mGapt.R1", s=5, alpha=0.5, color="cornflowerblue") #mRNA
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA1'], label = "T7-mGapt.R2", s=5, alpha=0.5, color="navy") #mRNA
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA2'], label = "T7-mGapt.R3", s=5, alpha=0.5, color="blue")#mRNA
plt.xlabel('Time (h)')
plt.ylabel('mGapt (uM)')
plt.legend()
R0.to_csv('MGapt_TXonly_PURE_simulation_old.csv', index=False) 
plt.savefig('MGapt_TXonly_PURE_simulation_old.png')


num_trajectories = 3
exp_data = []
for i in range(num_trajectories):
   df = pd.DataFrame()
   df["MGapt"]=truncated_data['mRNA'+str(i)]
   df["timepoints"] = truncated_data["time"]
   exp_data.append(df)
#exp_data

#Final parameters based on MGapt-TXwTL
k_rnapbF1 =10.194337
k_rnapbF2 = 5.270153
k_rnapbF3 = 10.425226
k_start=  7.819226
k_ntpbound=    2.451179
k_ntpadd=  97.083257
k_ntpdis=  1335.874684
k_term=  45.940369

# Choose the "true" parameters. NTP_add and NTP_deg
prior = {'k_forward' : ['gaussian', k_rnapbF1, 5, 'positive'],
        'k_forward1' : ['gaussian', k_rnapbF2, 2.5, 'positive'],
        'k_forward2' : ['gaussian', k_rnapbF3, 3, 'positive'], 
        'k_forward3' : ['gaussian', k_start, 4, 'positive'],
        'k_forward4' : ['gaussian', k_ntpbound, 1.5, 'positive'],
        'k_forward5' : ['gaussian', k_ntpadd, 30, 'positive'],
        'k_forward6' : ['gaussian', k_ntpdis, 300, 'positive'],
        'k_forward7' : ['gaussian', k_term, 15, 'positive'],
        }

sampler, pid = py_inference(Model = m, exp_data = exp_data, measurements = ['MGapt'],
                           time_column = ['timepoints'], 
                           params_to_estimate = ['k_forward', 'k_forward1', 'k_forward2',
                                                 'k_forward3', 'k_forward4','k_forward5','k_forward6','k_forward7',],
                           nwalkers = 30, nsteps = 500, init_seed = 'prior', prior = prior,
                           sim_type = 'deterministic', plot_show = True,debug=False, rtol=.01) #atol=.01



import pandas as pd
labels = ['k_forward', 'k_forward1', 'k_forward2','k_forward3', 
          'k_forward4', 'k_forward5','k_forward6',
         'k_forward7',]
samples_mcmc = pd.read_csv('mcmc_results.csv', names = labels,
                           engine = 'python')

import corner
import matplotlib.pyplot as plt
corner.corner(samples_mcmc, labels = labels, levels=(0.75,))
plt.savefig("MGapt_TXonly_corner_plot_March2024.svg")

#Get the mean ends
param_value=samples_mcmc.tail(500).mean()
print(samples_mcmc.tail(500).mean())

# Simulate model of mean
M_fit = Model(sbml_filename = filename+"_updated.xml")
for pi, pi_val in zip(['k_forward', 'k_forward1', 'k_forward2','k_forward3', 'k_forward4', 'k_forward5','k_forward6','k_forward7',], param_value):
        M_fit.set_parameter(pi, pi_val)
M_fit.set_species(initial_con)

R0 = py_simulate_model(Model = M_fit, timepoints = timepoints)
R0.to_csv('MGapt_TXonly_PURE_simulation_new.csv', index=False) 
plt.figure(figsize=(6, 8))
plt.title("mRNA Concentration Over Time")
plt.plot(R0['time']/3600, R0['MGapt'], label = "mRNA", color = 'm') #mRNA

#Experimental data set 1
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA0'], label = "T7-mGapt.R1", s=5, alpha=0.5, color="cornflowerblue") #mRNA
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA1'], label = "T7-mGapt.R2", s=5, alpha=0.5, color="navy") #mRNA
plt.scatter(truncated_data['time']/3600, truncated_data['mRNA2'], label = "T7-mGapt.R3", s=5, alpha=0.5, color="blue")#mRNA
plt.xlabel('Time (h)')
plt.ylabel('mGapt (uM)')
plt.legend()
plt.savefig('MGapt_TXonly_PURE_simulation_new.png')