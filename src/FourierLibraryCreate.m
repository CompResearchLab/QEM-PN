% generate POD modes (to run this script you need to make sure that main directory in in the path)
clear;
clc;
close all;
Input;Constants;
load("./Library/DataManagerPOD.mat") % loads the Trainer
% now that the Trainer is loaded compute the POD modes
DataM.FFTLib(NFML,FMM)
DataM.PSinvLib()
save("./Library/DataManagerPOD.mat","DataM",'-v7.3');


