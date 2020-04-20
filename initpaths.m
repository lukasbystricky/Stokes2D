%% Sets all the correct paths. 
% Assumes the git repository FastTools2D has been cloned. If not, it can be 
% found at https://github.com/lukasbystricky/FastTools2D.git

restoredefaultpath

close all
clear vars
clc

addpath(genpath('../FastTools2D/FMM/src/'))
addpath(genpath('../FastTools2D/PeriodicEwald/src'))
addpath(genpath('src'));
addpath(genpath('examples'));