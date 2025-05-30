UNS=1; % Number of random structures
NBS=1; % number of blocks (this code only runs with 1)
NS=1; % number of training structures
SingleG=SingleGauss(H,W,Shift,MinRadius,MaxRadius,Depth,n); % create gaussian block   
Blocks=cell(NBS,1); %define the number of blocks
Blocks{1}=SingleG; % set blocks