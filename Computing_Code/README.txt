# Repository for [The risk of SARS-CoV-2 outbreaks in low prevalence settings following the removal of travel restrictions](Published in Communications Medicine) #

Code author: Rahil Sachak-Patwa, University of Oxford
Last edit: 02/06/21
Software: MATLAB

System requirements: 
https://www.mathworks.com/support/requirements/matlab-system-requirements.html

Code tested on:
- Software Release: MATLAB_R2020a
- Operating system: macOS Big Sur Version 11.0.1
- Processor: 1.4 GHz Quad-Core Intel Core i5
- Memory: 16 GB 2133 MHx LPDD3

MATLAB installation link: https://www.mathworks.com/help/install/install-products.html
MATLAB Install time: ~2 hours

File descriptions
Data_IoM.mat: Matlab data file with data on historical and future projected vaccinations in the Isle of Man
Data_Israel.mat: Matlab data file with data on historical and future projected vaccinations in Israel
COR_ode: Ordinary differential equation function to compute the case outbreak risk (COR)
NOR_func.m: Function to compute the numerical outbreak risk (NOR)
outbreak_risk.m: Function to compute all four outbreak risk metrics, which calls NOR_func.m and COR_ode.m
fig_2_func.m: Function to plot the proportion of population vaccinated and compute the reproduction number over time, to plot Figure 2 in the manuscript 
fig_3_func.m: Function to plot the outbreak risk metrics, which calls outbreak_risk.m, to plot Figure 3 in the manuscript  
fig_settings.m: Function which defines the text formatting for the figures 
create_figures: Script which calls fig_2_func.m and fig_3_func.m which plots Figure 2 and Figure 3 in the manuscript  

Instructions:
- Extract all files from the zip file and add all files to a single folder.
- Add this folder to the MATLAB path.
- Run the MATLAB script create_figures.m which will produce Figures 2 and 3 in the main text. Note that the parameter values sim = 1e3 and P = 1 are set as default to increase the computational speed of the simulations, which will result in minor differences between the produced figures and the figures in the manuscript. 

Expected run time on a "normal" desktop computer: ~20 minutes.

This project is licensed under the terms of the MIT license.




