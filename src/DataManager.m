%TODO: allow for different size blocks
classdef DataManager <handle
    %This class is used for handeling the training
    
    properties
        % Training Data/POD data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Blocks
        Structure % function handle of the potential
        x % x data
        X % X data mesh
        y % y data mesh
        Y % Y data mesh
        n % Total Number of grid points per block
        nx % nodes in the x direction
        ny % nodes in the y direction
        dA
        BH % height of each block
        BW % width of each block
        rows
        cols
        CNT % stores how many structures have been trained for
        Configs % configure indicies 
        EEPOD
        a
        %Verification Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StructureTest
        asplit
        ModesI
        Vars
        TU
        Trows
        Tcols
        Tx
        Ty
        TX
        TY
        Tnx
        Tny
        DNSWF
        DNSEE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DNSTIME
        HamGTime
        HamSTime
        AssembleTimes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        dict % OFFDiagonal dictonary
    end

    methods
        % Blocks needs to be a cell!
        % n is the number of grid points per block
        function obj = DataManager(Blocks,rows,cols,n)
                % Get the total width and the total height 
                InitialB=Blocks{1}; % get the first block
                Height=InitialB.H; % get the height
                Width=InitialB.W; % get the Width
                obj.BH=Height;
                obj.BW=Width;
                obj.n=n;
                x=linspace(0,Width,n); % local x cordinates
                y=linspace(0,Height,n); % local y cordinates
                for i=1:cols-1
                    x=[x,x(2:n)+Width*i];
                end
                obj.nx=length(x);
                for i=1:rows-1
                    y=[y,y(2:n)+Height*i];
                end
                obj.ny=length(y);
                [X,Y]=meshgrid(x,y);
                obj.x=x;
                obj.y=y;
                dx=obj.x(2)-obj.x(1);
                dy=obj.y(2)-obj.y(1);
                obj.dA=dx*dy;
                obj.X=X;
                obj.Y=Y; 
                Indicies=zeros(rows*cols,1);
                obj.Blocks=Blocks;
                obj.Structure=cell(rows,cols);
                obj.rows=rows;
                obj.cols=cols;
                for i=1:rows
                    for j=1:cols
                        I=randi(length(Blocks));
                        obj.Structure{i,j}=Blocks{I};
                       
                        Indicies((i-1)*rows+j)=I;
                    end
                end
                obj.Configs=Indicies;
                obj.CNT=1;
        end
        % adjust
        function [] = TrainU(obj,rows,cols,UNRS,NRS)
                obj.rows=rows;
                obj.cols=cols;
                x=linspace(0,obj.BW,obj.n); % local x cordinates
                y=linspace(0,obj.BH,obj.n); % local y cordinates
                for i=1:cols-1
                    x=[x,x(2:obj.n)+obj.BW*i];
                end
                obj.nx=length(x);
                for i=1:rows-1
                    y=[y,y(2:obj.n)+obj.BH*i];
                end
                obj.ny=length(y);
                [X,Y]=meshgrid(x,y);
                obj.x=x;
                obj.y=y;
                obj.X=X;
                obj.Y=Y;
                
                obj.Structure=cell(rows,cols);
                for S=1:UNRS
                    for i=1:rows
                        for j=1:cols
                            I=randi(length(obj.Blocks));
                            obj.Structure{i,j}=obj.Blocks{I};
                        end
                    end
                    for RSI=1:NRS 
                        U=obj.generateU(false); % generate random sample variants of the structure
                        %savefig("./Library/"+string(S)+"_"+string(RSI))
                        obj.StoreUData(U);
                    end
                end
        end

        function []=StoreUData(obj,U)
             for i =1:obj.rows
                sy=(i-1)*obj.n-(i-2);
                ey=i*obj.n-(i-1);
                for j=1:obj.cols
                    sx=(j-1)*obj.n-(j-2);
                    ex=j*obj.n-(j-1);
                    CB=obj.Structure{i,j}; % current block
                    CU=U(sy:ey,sx:ex); % grab current potential
                     if isempty(CB.UData)
                            CB.UData=CU;
                     else
                            CB.UData(:,:,size(CB.UData,3)+1)=CU;
                     end
                end
             end
        end

        % add a function to generate the potential of the structure
        function U=generateU(obj,Flag)
            % Flag is to generate a contour if that is your choosing
            U=zeros(length(obj.x),length(obj.y)); % generate U;
            for i=1:obj.rows
                    y0=(i-1)*obj.BH;
                for j=1:obj.cols
                    x0=(j-1)*obj.BW;
                    CB=obj.Structure{i,j};
                    U=U+CB.genUR(obj.X,obj.Y,x0,y0); % apply superposition of potentials
                end
            end
            U=U-min(U,[],"all");
            if Flag==true
               obj.createUplot(obj.X,obj.Y,U,obj.cols,obj.rows);
            end
        end
        function []=createUplot(obj,X,Y,U,cols,rows)
             figure;
                sur=surf(X*1E9,Y*1E9,U);
                ax=gca;
                sur.EdgeAlpha=.1;
                x.ZLim(2)=ax.ZLim(2)*(1+.4);
                hold on
                con=contourf(X*1E9,Y*1E9,U,"ZLocation","zmax","FaceAlpha",.7,"EdgeAlpha",.7);
                %contourf(X*1E9,Y*1E9,U);
                MinZ=ax.ZLim(2);
                XM=ax.XLim(2);
                Xm=ax.XLim(1);
                YM=ax.YLim(2);
                Ym=ax.YLim(1);
                hold on
                plot3([Xm,Xm],[Ym,YM],[MinZ,MinZ],'k','linewidth',2);
                
                plot3([XM,XM],[Ym,YM],[MinZ,MinZ],'k','LineWidth',2);
                
                 for j=1:cols-1
                     plot3([j*obj.BW,j*obj.BW]*1E9,[Ym,YM],[MinZ,MinZ],'k','LineWidth',2)
                 end
                
                 plot3([Xm,XM],[Ym,Ym],[MinZ,MinZ],'k','linewidth',2);
               % plot([Xm,XM],[Ym,Ym],'k','linewidth',2);
                % 
                plot3([Xm,XM],[YM,YM],[MinZ,MinZ],'k','LineWidth',2);
                %plot([Xm,XM],[YM,YM],'k','LineWidth',2);
                 for j=1:rows-1
                     plot3([Xm,XM],[j*obj.BH,j*obj.BH]*1E9,[MinZ,MinZ],'k','LineWidth',2)
                 end
                % for j=1:rows-1
                %    plot([Xm,XM],[j*obj.BH,j*obj.BH]*1E9,'k','LineWidth',2)
                % end
                xlabel("X Position (nm)","FontSize",12);
                ylabel("Y Position (nm)","FontSize",12);
                zlabel("Energy (eV)","FontSize",12);

        end

        % This function mixes the current structure around
        function []=mix(obj)
            eq=true;
            while eq==true
                Indicies=zeros(obj.rows*obj.cols,1);
                for i=1:obj.rows
                    for j=1:obj.cols
                          I=randi(length(obj.Blocks));
                          Indicies((i-1)*obj.rows+j)=I;
                    end
                end
                % check to see if this structure has been used before 
                % if it has been used before we need to rerun the random
                % number generator.
                for i=1:size(obj.Configs,2)
                     if (isequal(Indicies,obj.Configs(:,i)))
                        eq=true;
                        break;
                     else
                         eq=false;
                     end
                end

            end
            obj.CNT=obj.CNT+1;
            obj.Configs(:,obj.CNT)=Indicies; % store the indicies
            for i =1:obj.rows
                for j=1:obj.cols
                    obj.Structure{i,j}=obj.Blocks{Indicies((i-1)*obj.rows+j)}; %rewrite the new structure
                end
            end
            %now that we have a unique structure we need to generate blocks
            disp(obj.Structure)
        end
        % Set the structure manually as a cell of blocks
        function []=setTest(obj,StructureTest)
            obj.StructureTest=cell(size(StructureTest,1),size(StructureTest,2)); % allocate memory
            obj.Vars=cell(size(StructureTest,1),size(StructureTest,2));
            obj.Trows=size(StructureTest,1);
            obj.Tcols=size(StructureTest,2);
            for i=1:size(StructureTest,1)
                for j=1:size(StructureTest,2)
                    CCell=StructureTest{i,j};
                    CB=CCell{1};
                    CV=CCell{2};
                    obj.StructureTest{i,j}=CB;
                    obj.Vars{i,j}=CV;
                end
            end
                Tx=linspace(0,obj.BW,obj.n); % local x cordinates
                Ty=linspace(0,obj.BH,obj.n); % local y cordinates
                for i=1:obj.Tcols-1
                    Tx=[Tx,Tx(2:obj.n)+obj.BW*i];
                end
                obj.Tnx=length(Tx);
                for i=1:obj.Trows-1
                    Ty=[Ty,Ty(2:obj.n)+obj.BH*i];
                end
                obj.Tny=length(Ty);
                [TX,TY]=meshgrid(Tx,Ty);
                obj.Tx=Tx;
                obj.Ty=Ty;
                obj.TX=TX;
                obj.TY=TY; 

        end
        function TU=generateUR(obj,rows,cols,Flag)
            obj.Trows=rows;
            obj.Tcols=cols;
            Tx=linspace(0,obj.BW,obj.n); % local x cordinates
            Ty=linspace(0,obj.BH,obj.n); % local y cordinates
            for i=1:obj.Tcols-1
                Tx=[Tx,Tx(2:obj.n)+obj.BW*i];
            end
            obj.Tnx=length(Tx);
            for i=1:obj.Trows-1
                Ty=[Ty,Ty(2:obj.n)+obj.BH*i];
            end
            obj.Tny=length(Ty);
            [TX,TY]=meshgrid(Tx,Ty);
            obj.Tx=Tx;
            obj.Ty=Ty;
            obj.TX=TX;
            obj.TY=TY; 
           
            obj.StructureTest=cell(rows,cols);
            obj.rows=rows;
            obj.cols=cols;
            for i=1:rows
                for j=1:cols
                    I=randi(length(obj.Blocks));
                    obj.StructureTest{i,j}=obj.Blocks{I};
                end
            end
            TU=zeros(length(obj.Tx),length(obj.Ty)); % generate U;
            for i=1:obj.Trows
                    y0=(i-1)*obj.BH;
                for j=1:obj.Tcols
                    x0=(j-1)*obj.BW;
                    CB=obj.StructureTest{i,j};
                    NewU=CB.genUR(obj.TX,obj.TY,x0,y0);
                    TU=TU+NewU; % apply superposition of potentials
                end
            end
            TU=TU-min(TU,[],'all');
            obj.TU=TU;

            if Flag==true
               obj.createUplot(obj.TX,obj.TY,TU,obj.Tcols,obj.Trows);
                
            end
        end


        function TU=genUTest(obj,Flag)
            TU=zeros(length(obj.Tx),length(obj.Ty)); % generate U;
            for i=1:obj.Trows
                    y0=(i-1)*obj.BH;
                for j=1:obj.Tcols
                    x0=(j-1)*obj.BW;
                    CB=obj.StructureTest{i,j};
                    CV=obj.Vars{i,j};
                    NewU=CB.genUT(obj.TX,obj.TY,x0,y0,CV);
                    TU=TU+NewU; % apply superposition of potentials
                end
            end
            TU=TU-min(TU,[],'all');
            obj.TU=TU;
             if Flag==true
                 obj.createUplot(obj.TX,obj.TY,TU,obj.Tcols,obj.Trows);
            end
        end

        function []=PlotPP(obj,Drop,NFM,PS,DoM,UMU)

            for i =1:PS:obj.rows
                sy=(i-1)*obj.n-(i-2);
                ey=i*obj.n-(i-1);
                for j=1:PS:obj.cols
                    sx=(j-1)*obj.n-(j-2);
                    ex=j*obj.n-(j-1);
                    CB=obj.Structure{i,j}; % current block
                    
                    CB.PlotPP(obj.TU(sy:ey,sx:ex),Drop,NFM,i,j,DoM,UMU);
                    
                   
                end
            end
        end
        
        % store the WF data corresponding to each block
        function []=StoreData(obj,WFS,U) % store WF and energy Data
           for i =1:obj.rows
                sy=(i-1)*obj.n-(i-2);
                ey=i*obj.n-(i-1);
                for j=1:obj.cols
                    sx=(j-1)*obj.n-(j-2);
                    ex=j*obj.n-(j-1);
                    CB=obj.Structure{i,j}; % current block
                    CU=U(sy:ey,sx:ex); % grab current potential
                     if isempty(CB.UData)
                            CB.UData=CU;
                     else
                            CB.UData(:,:,size(CB.UData,3)+1)=CU;
                     end
                    for s=1:size(WFS,2)
                        CW=reshape(WFS(:,s),obj.ny,obj.nx);
                        CWF=CW(sy:ey,sx:ex);
                        
                        if isempty(CB.WFData)
                            CB.WFData=CWF;
                        else
                            CB.WFData(:,:,size(CB.WFData,3)+1)=CWF;
                        end
                    end
                end
            end
        end
        function []=PlotAllEig(obj)
            figure;
            hold on
            for i=1:length(obj.Blocks)
                CB=obj.Blocks{i};
                semilogy(1:length(CB.EigenValues),CB.EigenValues,"LineWidth",2,"DisplayName",CB.Name);
            end
            set(gca,'yscale','log')
            hold off;
            legend();
            xlabel(" Mode Index ","FontSize",16)
            ylabel("Eigenvalue","FontSize",16)
            savefig("./Library/All_Eigs")
        end

       function []=FFTLib(obj,NFML,FMM)
  
            for i=1:length(obj.Blocks)
                CB=obj.Blocks{i};
                CB.UprepFFT(NFML,FMM)          
            end
            
       end
       function []=PSinvLib(obj)
             for i=1:length(obj.Blocks)
                CB=obj.Blocks{i};
                CB.PSinvLibCreate()
                CB.PSinvLibCreateMiddle()
                CB.PSinvLibCreateC()
                CB.PSinvLibCreateMiddleC()
            end
       end



        function []=PodLib(obj,LSize,Drop)
            
           
            for i=1:length(obj.Blocks)
                CB=obj.Blocks{i};
                CB.LSize=LSize;
                CB.ComputeModes();
                CB.ComputeUModes(Drop);
                CB.Kcomp();
                CB.Uprep();
                CB.DiagBcomp()
            end
            obj.PlotAllEig()
        end
          function []=ULib(obj,Drop)
            for i=1:length(obj.Blocks)
                CB=obj.Blocks{i};
                CB.ComputeUModes(Drop);
                CB.Uprep();
            end
            obj.PlotAllEig()
        end
        function []=MoveDataS(obj)
            % Move the data with the POD modes into the Structure Block
            obj.Structure=cell(obj.Trows,obj.Tcols); % store structure here
            for I=1:length(obj.Blocks)
                for i=1:obj.Trows
                    for j=1:obj.Tcols
                        if(strcmp(obj.Blocks{I}.Name,obj.StructureTest{i,j}.Name))
                            obj.Structure{i,j}=obj.Blocks{I};
                        end
                    end
                end
            end
        end
        function [LSE,NM,SingleNM,CompTime,TotalTime,LSEM]=compLSEU(obj,Drop,NFM,PM,MSize,nB,DoM,UMU)
            num=0;
            denom=0;
            NM=0;
            numM=0;
            denomM=0;
            TotalTime=0;
            for i =1:obj.Trows
                sy=(i-1)*obj.n-(i-2);
                ey=i*obj.n-(i-1);
                for j=1:obj.Tcols
                    sx=(j-1)*obj.n-(j-2);
                    ex=j*obj.n-(j-1);
                    CB=obj.Structure{i,j}; % current block
                    if(PM=="Classical")
                        [TNum,TDenom,nm,TNumM,TDenomM]=CB.ProjAEClassical(obj.TU(sy:ey,sx:ex),Drop,nB,DoM,UMU);
                        
                            tic
                            UM=CB.UCompClassical(obj.TU(sy:ey,sx:ex),MSize,Drop,DoM,UMU);
                            CompTime=toc;
                            TotalTime=TotalTime+CompTime;
                        
                    elseif(PM=="Fourier")
                        [TNum,TDenom,nm,TNumM,TDenomM]=CB.ProjAEFourier(obj.TU(sy:ey,sx:ex),NFM,nB);
                         
                            tic
                            UM=CB.UCompFourier(obj.TU(sy:ey,sx:ex),MSize,NFM);
                            CompTime=toc;
                            TotalTime=TotalTime+CompTime;

                        
                    elseif(PM=="PSinv")
                        [TNum,TDenom,nm,TNumM,TDenomM]=CB.ProjAEPSinv(obj.TU(sy:ey,sx:ex),Drop,nB,DoM,UMU);
                        
                            tic
                            UM=CB.UCompPSinv(obj.TU(sy:ey,sx:ex),MSize,Drop,DoM,UMU);
                            CompTime=toc;
                            TotalTime=CompTime+TotalTime;
                        
                    elseif(PM=="Normal")
                        
                            tic
                            UM=CB.UCompNormal(obj.TU(sy:ey,sx:ex),MSize);
                            LSE=0;
                            NM=0;
                            SingleNM=0;
                            CompTime=toc;
                            TotalTime=TotalTime+CompTime;
                            
                        
                    end
                    if(PM~="Normal")
                    if(i==1&&j==1)
                        SingleNM=nm;
                    end
                    NM=nm+NM;
                    numM=numM+TNumM;
                    denomM=denomM+TDenomM;
                    num=num+TNum;
                    denom=TDenom+denom; 
                    end
                end
            end
            if(PM~="Normal")
                LSE=sqrt(num/denom);
                LSEM=sqrt(numM/denomM);
            end
            
        end

        function H=genH(obj,MSize,N_r,Drop,NFM,PM,DoM,UMU)
            global PreOffDiag;
            NumBlocks=obj.Trows.*obj.Tcols;
            H=zeros(MSize*NumBlocks); % generate empty H matrix
            dx=obj.Tx(2)-obj.Tx(1);
            % load diagonals in the HMatrix
            for i =1:obj.Trows
                sy=(i-1)*obj.n-(i-2);
                ey=i*obj.n-(i-1);
                for j=1:obj.Tcols
                    sx=(j-1)*obj.n-(j-2);
                    ex=j*obj.n-(j-1);
                    CB=obj.Structure{i,j}; % current block
                    in=obj.GI(i,j); % get global index 
                    if(PM=="Classical")
                        Diag=(CB.UCompClassical(obj.TU(sy:ey,sx:ex),MSize,Drop,DoM,UMU)+CB.DT(1:MSize,1:MSize)+CB.DB1(1:MSize,1:MSize)+CB.DB2(1:MSize,1:MSize)*N_r/dx);  
                    elseif(PM=="Fourier")
                        Diag=(CB.UCompFourier(obj.TU(sy:ey,sx:ex),MSize,NFM)+CB.DT(1:MSize,1:MSize)+CB.DB1(1:MSize,1:MSize)+CB.DB2(1:MSize,1:MSize)*N_r/dx);  
                    elseif(PM=="PSinv")
                        Diag=(CB.UCompPSinv(obj.TU(sy:ey,sx:ex),MSize,Drop,DoM,UMU)+CB.DT(1:MSize,1:MSize)+CB.DB1(1:MSize,1:MSize)+CB.DB2(1:MSize,1:MSize)*N_r/dx); 
                    elseif(PM=="Normal")
                        Diag=(CB.UCompNormal(obj.TU(sy:ey,sx:ex),MSize)+CB.DT(1:MSize,1:MSize)+CB.DB1(1:MSize,1:MSize)+CB.DB2(1:MSize,1:MSize)*N_r/dx);  
                    end
                    H(1+MSize*(in-1):MSize*in,1+MSize*(in-1):MSize*in)=Diag;
                      
                end
            end

              % Compute diagonal entries off diagonal entries
            
             for i =1:obj.Trows
                for j=1:obj.Tcols
                    CB=obj.Structure{i,j}; % current block
                    inp=obj.GI(i,j); % get global index of present block
                    

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    TN=i+1; % Top neighbor
                    if(TN>obj.Trows)
                         TN=1;   
                    end
                    TNei=obj.Structure{TN,j};
                    inn=obj.GI(TN,j);
                    if(inn>inp)
                        if(PreOffDiag)
                            ODHC=obj.dict(string(CB.Name)+string(TNei.Name)+"T");
                            ODHC=ODHC{1}(1:MSize,1:MSize);
                        else
                        ODHC=obj.OFFDiag(CB,TNei,'T',MSize,N_r);
                        end
                        H(1+MSize*(inp-1):MSize*inp,1+MSize*(inn-1):MSize*inn)=ODHC;
                        H(1+MSize*(inn-1):MSize*inn,1+MSize*(inp-1):MSize*inp)=ODHC';
                    else
                        continue
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    BN=i-1; % Bottom Neighbor
                    if(BN==0)
                        BN=obj.Trows;
                    end
                    BNei=obj.Structure{BN,j};
                    inn=obj.GI(BN,j);
                    if(inn>inp)
                        if(PreOffDiag)
                            ODHC=obj.dict(string(CB.Name)+string(BNei.Name)+"B");
                            ODHC=ODHC{1}(1:MSize,1:MSize);
                        else
                        ODHC=obj.OFFDiag(CB,BNei,'B',MSize,N_r);
                        end
                        H(1+MSize*(inp-1):MSize*inp,1+MSize*(inn-1):MSize*inn)=ODHC;
                        H(1+MSize*(inn-1):MSize*inn,1+MSize*(inp-1):MSize*inp)=ODHC';
                    else
                        continue
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    RN=j+1; % Right Neighbor
                    if(RN>obj.Tcols)
                        RN=1;
                    end
                    RNei=obj.Structure{i,RN};
                    inn=obj.GI(i,RN);
                    if(inn>inp)
                        if(PreOffDiag)
                            ODHC=obj.dict(string(CB.Name)+string(RNei.Name)+"R");
                            ODHC=ODHC{1}(1:MSize,1:MSize);
                        else
                        ODHC=obj.OFFDiag(CB,RNei,'R',MSize,N_r);
                        end
                        H(1+MSize*(inp-1):MSize*inp,1+MSize*(inn-1):MSize*inn)= ODHC;
                        H(1+MSize*(inn-1):MSize*inn,1+MSize*(inp-1):MSize*inp)=ODHC';
                    else
                        continue
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    LN=j-1;
                    if(LN==0)
                        LN=obj.Tcols;
                    end
                    LNei=obj.Structure{i,LN};
                    inn=obj.GI(i,LN);
                    if(inn>inp)
                      if(PreOffDiag)
                      ODHC=obj.dict(string(CB.Name)+string(LNei.Name)+"L");
                      ODHC=ODHC{1}(1:MSize,1:MSize);
                      else
                      ODHC=obj.OFFDiag(CB,LNei,'L',MSize,N_r);
                      end
                      H(1+MSize*(inp-1):MSize*inp,1+MSize*(inn-1):MSize*inn)=ODHC;
                      H(1+MSize*(inn-1):MSize*inn,1+MSize*(inp-1):MSize*inp)=ODHC';
                    else
                        continue
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
             end
            H=(H'+H)/2;
        end
        function []=Solve(obj,MSize,N_r,CStates,Drop,NFM,PM,DoM,UMU,cutoffa)
             disp("Creating Hamiltonian")
             tic
             H=obj.genH(MSize,N_r,Drop,NFM,PM,DoM,UMU);
             time=toc
             obj.HamGTime=time;
             disp("Created Hamiltonian in: "+string(time) +"s")
             disp("Solving Hamiltonian")
             tic
             [a,EEPOD]=eigs(sparse(H),CStates,'sr');
             obj.EEPOD=diag(EEPOD);
             obj.a=a;
             amaxs=max(abs(a));
             obj.asplit=cell(obj.Trows,obj.Tcols,CStates);
             obj.ModesI=cell(obj.Trows,obj.Tcols,CStates);
               for i=1:obj.Trows 
                for j=1:obj.Tcols
                    inp=obj.GI(i,j); % get global index of present block
                    ac=obj.a(1+MSize*(inp-1):MSize*inp,:);
                    [as,MI]=sort(ac,1,'descend','ComparisonMethod','abs');
                    for s=1:CStates
                        I=find(abs(as(:,s)/amaxs(s))>cutoffa);
                        obj.asplit{i,j,s}=as(I,s);
                        obj.ModesI{i,j,s}=MI(I,s);
                    end
                   
                end
             end          
             time=toc
             obj.HamSTime=time;
             disp("Solving Hamiltonian in: "+string(time)+"s");
        end

        function DDictonary=ComputeOffDiag(obj,BP,BNei,LSize,N_r)
            Top=obj.OFFDiag(BP,BNei,'T',LSize,N_r);
            Bottom=obj.OFFDiag(BP,BNei,'B',LSize,N_r);
            Right=obj.OFFDiag(BP,BNei,'R',LSize,N_r);
            Left=obj.OFFDiag(BP,BNei,"L",LSize,N_r);
            keys=["T","B","R","L"];
            Matricies=cell(1,4);
            Matricies{1}=Top;
            Matricies{2}=Bottom;
            Matricies{3}=Right;
            Matricies{4}=Left;
            DDictonary=dictionary(keys,Matricies);
        end
        function []=PrecomputeOFFDiag(obj,LSize,N_r)
           NB=length(obj.Blocks); % get number of blocks in the library
           obj.dict=configureDictionary('string','cell');
          
           for i=1:NB
               BP=obj.Blocks{i};
                for j=1:NB
                    BNei=obj.Blocks{j};
                    DD=obj.ComputeOffDiag(BP,BNei,LSize,N_r);
                    obj.dict(string(BP.Name)+string(BNei.Name)+"T")=DD("T");
                    obj.dict(string(BP.Name)+string(BNei.Name)+"B")=DD("B");
                    obj.dict(string(BP.Name)+string(BNei.Name)+"R")=DD("R");
                    obj.dict(string(BP.Name)+string(BNei.Name)+"L")=DD("L");
                end
                
           end
        end
           


        function ODB=OFFDiag(obj,BP,BNei,Where,MSize,N_r)

             h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
             q_el = 1.60218e-19; % elementary charge, in coulombs
             Mass=9.10938215e-31; % Rest mass of particle in k
             C=h_bar^2/(2*Mass*q_el);
             x=linspace(0,BP.W,BP.n);
             y=linspace(0,BP.H,BP.n);
             dx=x(2)-x(1);
             dy=y(2)-y(1);
             mu=N_r/dx;
             ODB=zeros(MSize);
              for i=1:MSize
                switch(Where)
                    case 'T'
                        etai=BP.TM(:,:,i);
                        etai_dy=BP.TG(:,:,i);
                    case 'B'
                        etai=BP.BM(:,:,i);
                        etai_dy=BP.BG(:,:,i);
                    case 'R'
                        etai=BP.RM(:,:,i);
                        etai_dx=BP.RG(:,:,i);
                    case 'L'
                        etai=BP.LM(:,:,i);
                        etai_dx=BP.LG(:,:,i);
                end
                for j=1:MSize
                    switch(Where)
                        case 'T'
                            etaj=BNei.BM(:,:,i);
                            etaj_dy=BNei.BG(:,:,i);
                            F1=(etai_dy).*(etaj)-(etai).*(etaj_dy);
                            I=BP.xint(F1)/2;
                            F2=(etai).*(etaj);
                            I2=BP.xint(F2).*mu;
                            ODB(i,j)=C*(I-I2);
                        case 'B'
                            etaj=BNei.TM(:,:,j);
                            etaj_dy=BNei.TG(:,:,j);
                            F1=-((etai_dy).*(etaj)-(etai).*(etaj_dy));
                            I=BP.xint(F1)/2;
                            F2=(etai).*(etaj);
                            I2=BP.xint(F2)*mu;
                            ODB(i,j)=C*(I-I2);
                        case 'R'
                            etaj=BNei.LM(:,:,j);
                            etaj_dx=BNei.LG(:,:,j);
                            F1=(etai_dx).*(etaj)-(etai).*(etaj_dx);
                            I=BP.yint(F1)/2;
                            F2=(etai).*(etaj);
                            I2=BP.yint(F2)*mu;
                            ODB(i,j)=C*(I-I2);
                        case 'L'
                            etaj=BNei.RM(:,:,j);
                            etaj_dx=BNei.RG(:,:,j);
                            F1=-(etai_dx).*(etaj)-(etai).*(etaj_dx);
                            I=BP.yint(F1)/2;
                            F2=(etai).*(etaj);
                            I2=BP.yint(F2)*mu;
                            ODB(i,j)=C*(I-I2);
                        otherwise
                            exit("Something went Wrong")
                     end
                    
                end
              end
        end


        function in=GI(obj,i,j)
                in=(i-1)*obj.Tcols+j;
        end
      
   


    end
end

