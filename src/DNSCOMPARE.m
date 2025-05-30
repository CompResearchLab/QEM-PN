% DNS Comparision Runner (CODE ONLY WORKS FOR SQUARES)
clear;
clc;
close all;
Input; Constants; % run input file
load("./Library/DataManagerPOD.mat");
if(RANDt)
    UT=DataM.generateUR(Trows,Tcols,true);
else
    DataM.setTest(StructureTest);
    UT=DataM.genUTest(true); % generate potential
end
% figure;
% contourf(DataM.TX,DataM.TY,UT);
% savefig("./Results/U_Contour");
n=DataM.Tnx;
tic
B=FDTP(DataM.X,DataM.Y,n);
H=B+sparse(1:(n-1)*(n-1),1:(n-1)*(n-1),reshape(UT(1:end-1,1:end-1),[],1));
[rWFS,EE]= eigs(H,CStates,'sm'); % run eigenvalue solver
disp("DNS Solve and creation time is: ")
DataM.DNSTIME=toc;
DataM.DNSTIME
% put back the zeros
WFS=zeros(n,n,CStates);
for i=1:1:CStates
    WF=zeros(n,n);
    Z=reshape(rWFS(:,i),n-1,n-1);
    WF(1:end-1,:)=[Z,Z(:,1)];
    WF(end,:)=WF(1,:);
    WF=normalize(WF,DataM.dA)*WF; 
    WFS(:,:,i)=WF;    
end
DataM.DNSWF=WFS; % store the WFS in the prescribed class
DataM.DNSEE=diag(EE); 
disp("The DNS Energies are");
disp(DataM.DNSEE)
% note to load in the trainer in the POD library creation folder, one must
% also make sure the current file is added to the matlab PATH.
save("./Library/DataManagerPOD.mat","DataM",'-v7.3')


