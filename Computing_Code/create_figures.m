%% Create Figures 2 & 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Parmeters for Figures 2 & 3
S = 0.5;    % Standard deviation of R assuming normally distributed
mu = 7/5;   % Recovery rate (1/weeks)
del = 14;   % Delay in days between dose and susceptibility reduction 
t_end = 35;   % Number of weeks simulation is run for 
n = t_end*7;  % Number of days simulation is run for 
eta1 = 0.40; % Dose 1 efficacy parameter  
eta2 = 0.15; % Dose 2 efficacy parameter  

% Vector of countries and initial R values to run through to create figures
X_vec = [1,1,2,2];   % X = 1: IoM, X = 2: Israel
R_vec = [3,5,3,5]; % Initial R (reproduction number)

% Parmeters for Figure 3
M = 100;        % Threshold for outbreak (M)
sim = 1e3;      % Number of simulations in SOR, set as default to 1e3 for computational speed, increase to improve accuracy of SOR
P = 1;          % Number of samples values of R from normal distribution, set as default to 1 for computational speed, increase to sample more values of R

% Call function to plot Figure 2
fig_2_func(S,mu,del,n,eta1,eta2,X_vec,R_vec)

% Call function to plot Figure 3
tic
fig_3_func(S,mu,del,n,eta1,eta2,X_vec,R_vec,t_end,M,P,sim)
toc

