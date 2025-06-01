%INPUT FOR CODE
% Please note that because random numbers are involved you might need to
% adjust training and parameters based on the test case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************** Start of training data*************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh settings
n=200; % number of grid points per block in each direction
rows=2; % number of rows of the generated structures
cols=2; % number of columns of the generated structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train potentials (generate extra structures to train Potential)
UTrain=true; % specify true if you want extra training of the potentials
USizeR=5;% size of rows for training potentials
USizeC=5; % size of columns for training potentials
URS=60; % number of random samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=16E-9; % Height of each block (m)
W=16E-9; % Width of each block
Shift=0; % shift of the gaussian potential in each block
MaxRadius=6.5E-9; % max radius for each block
MinRadius=5E-9; % minimum radius for each block
Depth=0.8; %d epth of each guassian potential (eV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training settings
RS=50; % number of random samples for each structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SStates=15; % number of saved Modes 
SubSpace=100; % subspace for eigs for solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are the library settings (you need enough data to support these
% settings otherwise the code will not work)
USize=100; % Number of U Modes to train
LSize=100; % Size of POD Matrix Libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************End of Training Data**************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************Start Settings for POD-Galerkin evaluation and testing******
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CStates=15; % Number of collected states to test QEM model 
RANDt=true; % set true if you want a random comfiguration of elements
% If this value is set to false follow the template below to write your
% own.
Trows=4; % Total number of rows of elements
Tcols=4; % Total number of cols of elements

% If you want to generate a structure by hand follow the template below.
% The code was writen so that in the future other potentials besides
% gaussians can be used. However, this current version only supports
% gaussian. Hence you Must specify SingleG in each cell block. Then choose
% the shift of the gaussian potential in the x-coordinate, y-coordinate and
% the radius. All these values are in meters and the shift coordinate is in
% reference to the center of each element block.
% NOTE: UNCOMMENT AND EDIT IF YOU WHAT TO WRITE THE TEST STRUCTURE BY
% HAND!!!!
%StructureTest={{SingleG,[0,0,r(1)]},{SingleG,[0,0,r(5)]},{SingleG,[0,0,r(9)]},{SingleG,[0,0 r(13)]};...
%               {SingleG,[0,0,r(2)]},{SingleG,[0,0,r(6)]},{SingleG,[0,0,r(10)]},{SingleG,[0,0,r(14)]};...
%              {SingleG,[0,0,r(3)]},{SingleG,[0,0,r(7)]},{SingleG,[0,0,r(11)]},{SingleG,[0,0,r(15)]};
%                {SingleG,[0,0,r(4)]},{SingleG,[0,0,r(8)]},{SingleG,[0,0,r(12)]},{SingleG,[0,0,r(16)]} };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD Hamiltonian parameters
MSize=25; % Number of Modes for each element in Hamiltonian
N_r=.01; % pentalty factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%How do you want to specify the number of U-Modes (if you choose to use
%them)
DoM='C'; % Put D if you wan't to choose the number of Modes based on eigenvalue Drop, otherwise put 'C' for choose
UMU=15; % Number of U Modes used if using 'C'
Drop=1E-7; % Number of U modes used if selected 'D'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If DoM is set to 'D' and  if you are using POD U modes, plots showing the accuracy will appear when finished. The plot's x-axis is based
% off of order of magnitude drop from the first eigenvalue.
MinDrop=1; %x-axis starts at 1 order of magnitude drop
MaxDrop=9; %x-axis ends at 9 order of magnitude drop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine if you want to sweep across the penalty factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep across the penalty factor. This is done using a  logarithmically
% spaced scale.
N_RTest=false; % true if you want to preform sweep, elsewise false
N_Rmin=-2 % Starting decade of log scale
N_Rmax=1 % ending decade of log scale
N_Rn=30 % number of points to sweep across
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs if You want to sweep across Hamiltonian sizes
HTest=true; % Select true if you want to sweep across Hamiltonian sizes and plot the average error
HMin=20; % Minumum number of Modes used in Hamiltonian
HMax=40; % Maximum number of Modes used in Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************Potential Method**********************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program lets you choose the potential method as discussed in the
% paper. By projecting each the potential onto a know basis the
% U-Matrix for each element can be mostly precomputed before runtime.
% During runtime the only thing which needs to be computed are the basis
% coefficients.
PM="PSinv" % This variables selects the potential method. 
% Set it to "Normal" if you want to solve the potential matrix without
% any basis (real-space).
% Set it to "PSinv" or "Classical" if you want to project the potential onto POD U Modes.
% The coefficeints are then evaluated via projection or via LS fitting.
% Set it to "Fourier" if you want to project the potential onto fourier
% modes.
global NFML % if you select Fourier this global parameter is important
% The fourier library must be built before running QEMrunner.m.
NFML=15; % Maximum number of Fourier Modes in library
NFM=15; % Number of Fourier Modes used in simulation
FMM=30; % Number of POD Modes while building the library.
cutoffa=1E-3; % don't include a coefficients with a_i/max(a)<=cutoffa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSinv settings:
% If PSinv is set then you can adjust which points are included while
% using LS-fitting to find the U Mode Coefficients.
% The accuracy of the center of the potential is most import while using
% the QEM method. Hence, you can choose the to sample the middle points
% with greater accuracy than the edge.
global middle; 
middle=true; % Choose to sample the middle at a different ratio than points on the boundary.
global skip; % This parameter determines the point density of sampled points while computing the pseudo-inverse.
% choose an lower value to increase the accuracy.
skip=20; % Can only be 5,10,15,20,25 or 30

% Following two parameters only are taken into effect during training while
% creating the U-Library. If middle is true than you can select the number
% of boundary points and the spacing between sampled points.
global SB; SB=50; % Number of boundary points to skip (Only takes effect during training)
global Br; Br=20; %Boundary spacing resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting settings for the LSE plot regarding POD modes
contours=false % Plot contour plots
savefigsI=true %Save figures;
PlotScale=-.3; %Factor to Scale WFs to plot (-1, 0) Scales the y-axis by Depth*(1-PlotScale)
PlotStates=[6,9,11,15]
UsedModes=[1,7,10,0;1,9,17,0;1,9,20,0;1,10,18,0]; % used modes per state;
LSEM=20;% Number of Modes to include in LSE Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program will plot the potential in the basis you choose. It will plot
% the potential in each of the elements. This can generate a lot of plots
% if you are simulating a large structure. You can increase the plot
% incrementer to plot every PS elements.
PS=20; % Increment by Rows and columns by a greated number than 1 to reduce the number of Potential Error Plots
nB=20; % Since the center of the potential is the most important you can compute the Potential error of points nB:end-nB only.
%Set to 1 if you want the error for the entire potential.
   

