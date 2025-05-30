% DNS Training code runner
%TODO: CODE ONLY WORKS WITH dx=dy=dh (make it such that dx neq to dy)
clc;
clear;
close all;
Input;Constants;% read the input file
DataM=DataManager(Blocks,rows,cols,n);
x=DataM.x; % get x cordintes
y=DataM.y; % get y cordinates
X=DataM.X; % get X mesh
Y=DataM.Y; % get Y mesh
n=DataM.nx; % get the number of mesh points
B=FDTP(X,Y,n);

CNTER=1;
for structI=1:NS % loop through NS structures
   disp("Runing Structure "+string(structI)+" out of "+string(NS))
   if(structI~=1)
        DataM.mix(); % change the structure
        disp(DataM.Structure) % write out the structure to the terminal
   end  
   
    for RI=1:RS % loop through RS random iterations of the structure
            U=DataM.generateU(false); % get a new U
            %savefig("./Library/U_Contour"+string(CNTER))
            CNTER=CNTER+1;
            H=B+sparse(1:(n-1)*(n-1),1:(n-1)*(n-1),reshape(U(1:end-1,1:end-1),[],1));
            [rWFS,EE]= eigs(H,SStates,'sm'); % run eigenvalue solver
            % put back the zeros
            WFS=zeros(n*n,SStates);
            for i=1:1:SStates
                WF=zeros(n,n);
                Z=reshape(rWFS(:,i),n-1,n-1);
                WF(1:end-1,:)=[Z,Z(:,1)];
                WF(end,:)=WF(1,:);
                WF=normalize(WF,DataM.dA)*WF; 
                WFS(:,i)=reshape(WF,n*n,1);    
            end
            DataM.StoreData(WFS,U); % store the WFS in the prescribed class
           
            disp("Running variation "+ string(RI)+ " out of "+string(RS))
     end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is extra training if one further needs to train the Umodes
if(UTrain==true)
    DataM.TrainU(USizeR,USizeC,UNS,URS);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note to load in the trainer in the POD library creation folder, one must
% also make sure the current file is added to the matlab PATH.
save("./Library/DataManagerPOD","DataM",'-v7.3')


