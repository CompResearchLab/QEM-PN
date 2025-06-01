% generate POD modes (to run this script you need to make sure that main directory in in the path)
clear;
clc;
close all;
Input;
load("./Library/DataManagerPOD.mat") % loads the Trainer
DataM.MoveDataS();
DataM.Solve(MSize,N_r,CStates,Drop,NFM,PM,DoM,UMU,cutoffa);
%%
Post=PostProcess(DataM);

Post.PlotStates(PlotStates,UsedModes,CStates,MSize,N_Rmin,N_Rmax,N_Rn,N_RTest,savefigsI,contours,HTest,HMin,HMax,N_r,MaxDrop,MinDrop,PlotScale,LSEM,Drop,PM,NFM,PS,nB,UMU,DoM,cutoffa);