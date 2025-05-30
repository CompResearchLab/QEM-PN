classdef (Abstract) Block <handle
    % This is the block class

    properties(Access=public)
        Name % name of block must be a string
        W {mustBeNumeric} % width
        H {mustBeNumeric} % height
        n % Total Number of grid points in each direction for the block
        x
        y
        dx
        dy
        dA
        WFData % Wave function data;
        Modes % store the modes in a library
        ModesMult % Stores Mode in form to be multiplied
        UData
        UModes
        EigenValuesU
        UModesInt
        NumUM
        FX % ddeta/dx
        FY% deta/dy
        EigenValues
        LSize % size of POD Matricies
        DT % diagonal kinetic energy matrix
        DB1 % diagonal boundary matrix
        DB2
        DUprep
        DUprepFourier;
      
        % Inteface of Modes
        TM
        BM
        LM
        RM
        % Interface of gradients
        TG
        BG
        LG
        RG
        LSE
       % PSlib
       PSinv5;   PSinv5C;
       PSinv10;  PSinv10C;
       PSinv15;  PSinv15C;
       PSinv20;  PSinv20C;
       PSinv25;  PSinv25C;
       PSinv30;  PSinv30C;
       PSinvM5;  PSinvM5C;
       PSinvM10; PSinvM10C;
       PSinvM15; PSinvM15C;
       PSinvM20; PSinvM20C;
       PSinvM25; PSinvM25C;
       PSinvM30; PSinvM30C;
    end

    methods(Abstract)
       U=genUR(obj,X,Y,x0,y0)
       U=genUT(obj,X,Y,x0,y0,Vars)
    end
    methods
        % get top of Matrix
        function T=gT(obj,M)
            T=M(end,:);
        end
        % get bottom of Matrix
        function B=gB(obj,M)
            B=M(1,:); 
        end
        % get Right
        function R=gR(obj,M)
            R=M(:,end);
        end
        function L=gL(obj,M)
            L=M(:,1);
        end
        function obj=Block()
                obj.Name=class(obj); % get the name of the block and put it into
                % the Name property
        end
        function I=xint(obj,f)
            % integrate horizontal direction
            I=obj.dx.*((f(1)+f(end))./2+sum(f(2:end-1)));
        end
        function I=yint(obj,f)
           % integrate vertical direction
           I=obj.dy.*((f(1)+f(end))./2+sum(f(2:end-1)));
        end
        function I=BINT(obj,f)
           
            I=Int2D(f,obj.dA); % compute integral
        end
        function []=Kcomp(obj)
            %SnapShots=size(obj.Modes,3);
            disp("Computing DT for "+ obj.Name)
             % generate the kinetic energy operator in finite difference form
             h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
             q_el = 1.60218e-19; % elementary charge, in coulombs
             Mass=9.10938215e-31; % Rest mass of particle in kg
             T=zeros(obj.LSize);
         
            for i=1:obj.LSize
                for j=i:obj.LSize
                    T(i,j)=obj.BINT(h_bar^2./(2*Mass*q_el).*(obj.FX(:,:,i).*obj.FX(:,:,j)+obj.FY(:,:,i).*obj.FY(:,:,j)));
                end
            end
            % exploite symmetry
            for i=1:obj.LSize
                T(i,i)=T(i,i)/2;
            end
            T=T'+T;
            obj.DT=T;
        end


    function []=Uprep(obj)
       obj.DUprep=zeros(obj.LSize,obj.LSize,obj.NumUM);
       for umi=1:obj.NumUM
         umode=obj.UModes(:,:,umi); % grab the current mode
         for i=1:obj.LSize
            etai=obj.Modes(:,:,i); % WF modes
            for j=i:obj.LSize
                etaj=obj.Modes(:,:,j); % WF modes
                obj.DUprep(i,j,umi)=obj.BINT(etai.*etaj.*umode);
            end
         end
        % exploite symmetry
        for i=1:obj.LSize
            obj.DUprep(i,i,umi)=obj.DUprep(i,i,umi)./2;
        end
        obj.DUprep(:,:,umi)=obj.DUprep(:,:,umi)'+obj.DUprep(:,:,umi);
        obj.UModesInt(:,:,umi)=obj.UModes(:,:,umi).*obj.dA;
        obj.UModesInt(1,:,umi)=obj.UModesInt(1,:,umi)/2;
        obj.UModesInt(end,:,umi)=obj.UModesInt(end,:,umi)/2;
        obj.UModesInt(:,1,umi)=obj.UModesInt(:,1,umi)/2;
        obj.UModesInt(:,end,umi)=obj.UModesInt(:,end,umi)/2;
        % obj.UModesInt(1,1,umi)=obj.UModesInt(1,1,umi)/2;
        % obj.UModesInt(end,end,umi)=obj.UModesInt(end,end,umi)/2;
        % obj.UModesInt(1,end,umi)=obj.UModesInt(1,end,umi)/2;
        % obj.UModesInt(end,1,umi)=obj.UModesInt(end,1,umi)/2;
       end
    end


    function []=PSinvLibCreateMiddle(obj)
        global SB;global Br; % Boundary skip
        
        skip=[5,10,15,20,25,30];
        Drop=1:9;
        obj.PSinvM5=cell(length(Drop));
        obj.PSinvM10=cell(length(Drop));
        obj.PSinvM15=cell(length(Drop));
        obj.PSinvM20=cell(length(Drop));
        obj.PSinvM25=cell(length(Drop));
        obj.PSinvM30=cell(length(Drop));
        for i=Drop
            NM=find(obj.EigenValuesU<10^(-i),1,'first');
            for s=skip
                switch (s) 
                    case 5
                        obj.PSinvM5{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);

                    case 10
                        obj.PSinvM10{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 15
                        obj.PSinvM15{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 20
                        obj.PSinvM20{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 25
                        obj.PSinvM25{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 30
                        obj.PSinvM30{i}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                end
            end
        end
    end
    function []=PSinvLibCreateMiddleC(obj)
        global SB;global Br; % Boundary skip
        OrderDropMaxNumberOfModes=9; % order drop to get maximum number of modes
        NumMMax=find(obj.EigenValuesU<10^(-OrderDropMaxNumberOfModes),1,'first');% Maxnumber of modes
        skip=[5,10,15,20,25,30];
        
        obj.PSinvM5C=cell(1,NumMMax);
        obj.PSinvM10C=cell(1,NumMMax);
        obj.PSinvM15C=cell(1,NumMMax);
        obj.PSinvM20C=cell(1,NumMMax);
        obj.PSinvM25=cell(1,NumMMax);
        obj.PSinvM30C=cell(1,NumMMax);
        for NM=1:NumMMax
            
            for s=skip
                switch (s) 
                    case 5
                        obj.PSinvM5C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);

                    case 10
                        obj.PSinvM10C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 15
                        obj.PSinvM15C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 20
                        obj.PSinvM20C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 25
                        obj.PSinvM25C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                    case 30
                        obj.PSinvM30C{NM}=pinv([reshape(obj.UModes(SB:s:end-SB,SB:s:end-SB,1:NM),[],NM);... % boundary middle
                           reshape(obj.UModes(1:Br:end,1:Br:SB-1,1:NM),[],NM);... % boundary first set of columns
                           reshape(obj.UModes(1:Br:end,end-SB+1:Br:end,1:NM),[],NM);... % boundary second set of columns
                           reshape(obj.UModes(1:Br:SB-1,SB:Br:end-SB,1:NM),[],NM);
                           reshape(obj.UModes(end-SB+1:Br:end,SB:Br:end-SB,1:NM),[],NM);]);
                end
            end
        end
    end


    function []=PSinvLibCreateC(obj)
        skip=[5,10,15,20,25,30];
        OrderDropMaxNumberOfModes=9; % order drop to get maximum number of modes
        NumMMax=find(obj.EigenValuesU<10^(-OrderDropMaxNumberOfModes),1,'first');% Maxnumber of modes
        obj.PSinv5C=cell(1,NumMMax);
        obj.PSinv10C=cell(1,NumMMax);
        obj.PSinv15C=cell(1,NumMMax);
        obj.PSinv20C=cell(1,NumMMax);
        obj.PSinv25C=cell(1,NumMMax);
        obj.PSinv30C=cell(1,NumMMax);
        for NM=1:NumMMax
            
            for s=skip
                switch (s) 
                    case 5
                        obj.PSinv5{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 10
                        obj.PSinv10{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 15
                        obj.PSinv15{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 20
                        obj.PSinv20{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 25
                        obj.PSinv25{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 30
                        obj.PSinv30{NM}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                end
            end
        end
    end

    function []=PSinvLibCreate(obj)
        skip=[5,10,15,20,25,30];
        Drop=1:9;
        obj.PSinv5=cell(length(Drop));
        obj.PSinv10=cell(length(Drop));
        obj.PSinv15=cell(length(Drop));
        obj.PSinv20=cell(length(Drop));
        obj.PSinv25=cell(length(Drop));
        obj.PSinv30=cell(length(Drop));
        for i=Drop
            NM=find(obj.EigenValuesU<10^(-i),1,'first');
            for s=skip
                switch (s) 
                    case 5
                        obj.PSinv5{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 10
                        obj.PSinv10{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 15
                        obj.PSinv15{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 20
                        obj.PSinv20{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 25
                        obj.PSinv25{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                    case 30
                        obj.PSinv30{i}=pinv(reshape(obj.UModes(1:s:end,1:s:end,1:NM),[],NM));
                end
            end
        end
    end
    function DU=UCompNormal(obj,U,MSize)
        DU=zeros(MSize);
        for i=1:MSize
            etai=obj.Modes(:,:,i); % WF modes
            for j=i:MSize
                etaj=obj.Modes(:,:,j); % WF modes
                DU(i,j)=obj.BINT(etai.*etaj.*U);
            end
        end
        % exploite symmetry
        for i=1:MSize
            DU(i,i)=DU(i,i)./2;
        end
        DU=DU'+DU;
    end
function DU=UCompPSinv(obj,U,MSize,Drop,DoM,UMU)
       global skip; global middle; global SB; global Br;
       % grab weights
       DU=zeros(MSize);
       if(DoM=='D')
       DropI=-1*log10(Drop);
       
       if(middle==false)
        switch skip
           case 5
                M=obj.PSinv5{DropI};
           case 10
               M=obj.PSinv10{DropI};
           case 15
               M=obj.PSinv15{DropI};
           case 20
               M=obj.PSinv20{DropI};
           case 25
               M=obj.PSinv25{DropI};
           case 30
               M=obj.PSinv30{DropI};
        end
       else
        switch skip
           case 5
                M=obj.PSinvM5{DropI};
           case 10
               M=obj.PSinvM10{DropI};
           case 15
               M=obj.PSinvM15{DropI};
           case 20
               M=obj.PSinvM20{DropI};
           case 25
               M=obj.PSinvM25{DropI};
           case 30
               M=obj.PSinvM30{DropI};
        end
        end
       
       
       if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
           C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);         
       end
    NM=find(obj.EigenValuesU<Drop,1,'first');
    for i=1:NM
        DU=DU+C(i)*obj.DUprep(1:MSize,1:MSize,i);
    end
       else
            if(middle==false)
                switch skip
                case 5
                    M=obj.PSinv5C{UMU};
                case 10
                    M=obj.PSinv10C{UMU};
                case 15
                    M=obj.PSinv15C{UMU};
                case 20
                    M=obj.PSinv20C{UMU};
                case 25
                    M=obj.PSinv25C{UMU};
                case 30
                    M=obj.PSinv30C{UMU};
                end
       else
        switch skip
           case 5
                M=obj.PSinvM5C{UMU};
           case 10
               M=obj.PSinvM10C{UMU};
           case 15
               M=obj.PSinvM15C{UMU};
           case 20
               M=obj.PSinvM20C{UMU};
           case 25
               M=obj.PSinvM25C{UMU};
           case 30
               M=obj.PSinvM30C{UMU};
        end
        end
       
       
       if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
           C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);         
       end
    
    for i=1:UMU
        DU=DU+C(i)*obj.DUprep(1:MSize,1:MSize,i);
    end
  end
end



function []=UprepFFT(obj,NFML,FMM)
      obj.DUprepFourier=zeros(obj.LSize,obj.LSize,NFML,NFML); 
      x=linspace(0,obj.W,obj.n);
      y=linspace(0,obj.H,obj.n);
      [X,Y]=meshgrid(x,y);
      n=obj.n; ;W=obj.W; H=obj.H;
     if(mod(n,2)==0) % even split
        dk=[-n/2:1:1:n/2-1]; c=n/2+1; % <-- using n
     else %odd split
        dk=[-(n-1)/2:1:(n-1)/2]; c=n/2; % <-- using n
     end
     if(mod(NFML,2)==0) % even split
        dkN=[-NFML/2:1:1:NFML/2-1]; % <-- using n
     else %odd split
        dkN=[-(NFML-1)/2:1:(NFML-1)/2]; % <-- using n
     end
      [Xn,Yn]=meshgrid([0:n-1],[0:n-1]); % mesh for FFT
      wm=exp(2*pi*1i/n); % exponentials for FFT
      wn=exp(2*pi*1i/n); %exponentials for FFT
      for jk=(1:NFML)
        for ik=(1:NFML)
            e=wm.^(Yn.*dkN(jk)).*wn.^(dkN(ik).*Xn);
            for i=1:FMM
                etai=obj.Modes(:,:,i); % WF modes
                for j=i:FMM
                    etaj=obj.Modes(:,:,j); % WF modes
                    obj.DUprepFourier(i,j,jk,ik)=obj.BINT(etai.*e.*etaj);
                end
            
            end
            for z=1:FMM
            obj.DUprepFourier(z,z,jk,ik)=obj.DUprepFourier(z,z,jk,ik)/2;
            end
            obj.DUprepFourier(:,:,jk,ik)=obj.DUprepFourier(:,:,jk,ik).'+obj.DUprepFourier(:,:,jk,ik);
         end
      end 
end
function DU=UCompClassical(obj,U,MSize,Drop,DoM,UMU)
            % grab weights
            DU=zeros(MSize);
           % UP=zeros(obj.n,obj.n);
           if(DoM=="D")
            NM=find(obj.EigenValuesU<Drop,1,'first');
            for i=1:NM
                C=obj.BINT(U.*obj.UModes(:,:,i));
                %UP=UP+C*obj.UModes(:,:,i);
                DU=DU+C*obj.DUprep(1:MSize,1:MSize,i); 
            end   
           else
               for i=1:UMU
                C=obj.BINT(U.*obj.UModes(:,:,i));
                DU=DU+C*obj.DUprep(1:MSize,1:MSize,i); 
                end 

           end
    end

    
  function DU=UCompFourier(obj,U,MSize,NFM)
 
        global NFML; % number of modes in library
        n=obj.n;
        if(mod(n,2)==0) % even split
            dk=[-n/2:1:1:n/2-1]; c=n/2+1; % <-- using n
        else %odd split
            dk=[-(n-1)/2:1:(n-1)/2]; c=ceil(n/2); % <-- using n
        end
        if(mod(NFM,2)==0) % even split
            dkN=[-NFM/2:1:1:NFM/2-1]; % <-- using n
        else %odd split
            dkN=[-(NFM-1)/2:1:(NFM-1)/2]; % <-- using n
        end
        if(mod(NFML,2)==0) % even split
            dkN=[-NFML/2:1:1:NFML/2-1];cl=NFML/2+1; % <-- using n
        else %odd split
            dkN=[-(NFML-1)/2:1:(NFML-1)/2];cl=ceil(NFML/2); % <-- using n
        end

        % grab weights
       Coefs=1/obj.n^2*fftshift(fft2(U));
       
       DU=zeros(MSize);
       for j=1:NFM
            for k=1:NFM
                DU=DU+Coefs(c+dkN(j),c+dkN(k)).*obj.DUprepFourier(1:MSize,1:MSize,dkN(j)+cl,dkN(k)+cl);
            end
       end
       
       DU=real(DU);       
  end
  function [num,denom,NM,numM,denomM]=ProjAEPSinv(obj,U,Drop,nB,DoM,UMU)
        global skip; global middle; global SB; global Br;
      
        
        UP=zeros(obj.n,obj.n);
       if(DoM=='D')
       DropI=-1*log10(Drop);
       NM=find(obj.EigenValuesU<Drop,1,'first');
       if(middle==false)
        switch skip
           case 5
                M=obj.PSinv5{DropI};
           case 10
               M=obj.PSinv10{DropI};
           case 15
               M=obj.PSinv15{DropI};
           case 20
               M=obj.PSinv20{DropI};
           case 25
               M=obj.PSinv25{DropI};
           case 30
               M=obj.PSinv30{DropI};
        end
       else
        switch skip
           case 5
                M=obj.PSinvM5{DropI};
           case 10
               M=obj.PSinvM10{DropI};
           case 15
               M=obj.PSinvM15{DropI};
           case 20
               M=obj.PSinvM20{DropI};
           case 25
               M=obj.PSinvM25{DropI};
           case 30
               M=obj.PSinvM30{DropI};
        end
        end
       
        
        
        if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
             C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);
       end
        for i=1:NM
                UP=UP+C(i)*obj.UModes(:,:,i);
        end
       else

       NM=UMU;
       if(middle==false)
        switch skip
           case 5
                M=obj.PSinv5C{UMU};
           case 10
               M=obj.PSinv10C{UMU};
           case 15
               M=obj.PSinv15C{UMU};
           case 20
               M=obj.PSinv20C{UMU};
           case 25
               M=obj.PSinv25C{UMU};
           case 30
               M=obj.PSinv30C{UMU};
        end
       else
        switch skip
           case 5
                M=obj.PSinvM5C{UMU};
           case 10
               M=obj.PSinvM10C{UMU};
           case 15
               M=obj.PSinvM15C{UMU};
           case 20
               M=obj.PSinvM20C{UMU};
           case 25
               M=obj.PSinvM25C{UMU};
           case 30
               M=obj.PSinvM30C{UMU};
        end
        end
       
        
        
        if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
             C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);
       end
        for i=1:NM
                UP=UP+C(i)*obj.UModes(:,:,i);
        end
       end
        x=linspace(0,obj.W,obj.n);
        y=linspace(0,obj.H,obj.n);
        [X,Y]=meshgrid(x,y);
        num=obj.BINT((UP-U).^2);
        denom=obj.BINT((U).^2);
        M=(UP-U).^2;
        numM=obj.BINT(M(nB:end-nB,nB:end-nB));
        denomM=obj.BINT(U(nB:end-nB,nB:end-nB).^2);

    end

    function [num,denom,NM,numM,denomM]=ProjAEClassical(obj,U,Drop,nB,DoM,UMU)
        if(DoM=="D")
        NM=find(obj.EigenValuesU<Drop,1,'first');
        UP=zeros(obj.n,obj.n);
        for i=1:NM
                C=sum(U.*obj.UModesInt(:,:,i),'all');
                UP=UP+C*obj.UModes(:,:,i);
        end 
        else
        NM=UMU;
        UP=zeros(obj.n,obj.n);
        for i=1:UMU
                C=sum(U.*obj.UModesInt(:,:,i),'all');
                UP=UP+C*obj.UModes(:,:,i);
        end 
        end
        x=linspace(0,obj.W,obj.n);
        y=linspace(0,obj.H,obj.n);
        [X,Y]=meshgrid(x,y);
        num=obj.BINT((UP-U).^2);
        denom=obj.BINT((U).^2);
        M=(UP-U).^2;
        numM=obj.BINT(M(nB:end-nB,nB:end-nB));
        denomM=obj.BINT(U(nB:end-nB,nB:end-nB).^2);


  end
    function []=PlotPP(obj,U,Drop,NFM,row,col,DoM,UMU)
        global skip; global middle; global SB; global Br;
        NM=find(obj.EigenValuesU<Drop,1,'first');
        UP=zeros(obj.n,obj.n);
        if(DoM=="D")
        for i=1:NM
                C=obj.BINT(U.*obj.UModes(:,:,i));
                
                UP=UP+C*obj.UModes(:,:,i);
        end
        else
        for i=1:UMU
                C=obj.BINT(U.*obj.UModes(:,:,i));
                
                UP=UP+C*obj.UModes(:,:,i);
        end
        end
        plt=figure; hold on;
        x=linspace(0,obj.W,obj.n);
        plot(x*1E9,U(floor(obj.n/2),:),'k','LineWidth',2,'DisplayName'," True Potential")
        plot(x*1E9,UP(floor(obj.n/2),:),'r:','LineWidth',3,'DisplayName'," Classical")
        xlabel("X Position (nm)",'FontSize',16)
        ylabel("Potential (eV)","FontSize",16)
        % Now Look at PSinv Method
        UP=zeros(obj.n,obj.n);
        if(DoM=="D")
        DropI=-1*log10(Drop);
       
       if(middle==false)
        switch skip
           case 5
                M=obj.PSinv5{DropI};
           case 10
               M=obj.PSinv10{DropI};
           case 15
               M=obj.PSinv15{DropI};
           case 20
               M=obj.PSinv20{DropI};
           case 25
               M=obj.PSinv25{DropI};
           case 30
               M=obj.PSinv30{DropI};
        end
       else
        switch skip
           case 5
                M=obj.PSinvM5{DropI};
           case 10
               M=obj.PSinvM10{DropI};
           case 15
               M=obj.PSinvM15{DropI};
           case 20
               M=obj.PSinvM20{DropI};
           case 25
               M=obj.PSinvM25{DropI};
           case 30
               M=obj.PSinvM30{DropI};
        end
        end
       
        
        
        if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
             C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);
       end
           
        for i=1:NM
                UP=UP+C(i)*obj.UModes(:,:,i);
        end
        else
             NM=UMU;
       if(middle==false)
        switch skip
           case 5
                M=obj.PSinv5C{UMU};
           case 10
               M=obj.PSinv10C{UMU};
           case 15
               M=obj.PSinv15C{UMU};
           case 20
               M=obj.PSinv20C{UMU};
           case 25
               M=obj.PSinv25C{UMU};
           case 30
               M=obj.PSinv30C{UMU};
        end
       else
        switch skip
           case 5
                M=obj.PSinvM5C{UMU};
           case 10
               M=obj.PSinvM10C{UMU};
           case 15
               M=obj.PSinvM15C{UMU};
           case 20
               M=obj.PSinvM20C{UMU};
           case 25
               M=obj.PSinvM25C{UMU};
           case 30
               M=obj.PSinvM30C{UMU};
        end
        end
       
        
        
        if(middle==false)
            C=M*reshape(U(1:skip:end,1:skip:end),[],1);
       else
             C=M*([reshape(U(SB:skip:end-SB,SB:skip:end-SB),[],1);... % boundary middle
                           reshape(U(1:Br:end,1:Br:SB-1),[],1);... % boundary first set of columns
                           reshape(U(1:Br:end,end-SB+1:Br:end),[],1);... % boundary second set of columns
                           reshape(U(1:Br:SB-1,SB:Br:end-SB),[],1);
                           reshape(U(end-SB+1:Br:end,SB:Br:end-SB),[],1);]);
       end
        for i=1:NM
                UP=UP+C(i)*obj.UModes(:,:,i);
        end


        end
        plot(x*1E9,UP(floor(obj.n/2),:),'b--','LineWidth',3,'DisplayName'," Least Squares")
        % Now look at Fourier Method
        n=obj.n;
        if(mod(n,2)==0) % even split
            dk=[-n/2:1:1:n/2-1]; c=n/2+1; % <-- using n
        else %odd split
            dk=[-(n-1)/2:1:(n-1)/2]; c=n/2; % <-- using n
        end
        if(mod(NFM,2)==0) % even split
            dkN=[-NFM/2:1:1:NFM/2-1]; % <-- using n
        else %odd split
            dkN=[-(NFM-1)/2:1:(NFM-1)/2]; % <-- using n
        end
        fhat=fftshift(fft2(U));
        wm=exp(2*pi*1i/n);
        wn=exp(2*pi*1i/n);
        UP=zeros(n);
        [Xn,Yn]=meshgrid([0:n-1],[0:n-1]);
        for j=1:NFM
            for k=1:NFM
                 UP=UP+fhat(c+dkN(j),c+dkN(k)).*wm.^(Yn.*dkN(j)).*wn.^(dkN(k).*Xn);
            end
        end
        UP=real(UP*1/n.^2);
        plot(x*1E9,UP(floor(obj.n/2),:),'g-.','LineWidth',2,'DisplayName'," Fourier")
        legend("FontSize",16)
        savefig(plt,"./Results/ProfilePlotPotentialError_"+string(row)+"_"+string(col));
    end

    function [num,denom,NM,numM,denomM]=ProjAEFourier(obj,U,NFM,nB)
        n=obj.n;
        if(mod(n,2)==0) % even split
            dk=[-n/2:1:1:n/2-1]; c=n/2+1; % <-- using n
        else %odd split
            dk=[-(n-1)/2:1:(n-1)/2]; c=n/2; % <-- using n
        end
        if(mod(NFM,2)==0) % even split
            dkN=[-NFM/2:1:1:NFM/2-1]; % <-- using n
        else %odd split
            dkN=[-(NFM-1)/2:1:(NFM-1)/2]; % <-- using n
        end
        fhat=fftshift(fft2(U));
        wm=exp(2*pi*1i/n);
        wn=exp(2*pi*1i/n);
        UP=zeros(n);
        [Xn,Yn]=meshgrid([0:n-1],[0:n-1]);
        for j=1:NFM
            for k=1:NFM
                 UP=UP+fhat(c+dkN(j),c+dkN(k)).*wm.^(Yn.*dkN(j)).*wn.^(dkN(k).*Xn);
            end
        end
        UP=real(UP*1/n.^2);
        NM=NFM.*NFM;
        x=linspace(0,obj.W,obj.n);
        y=linspace(0,obj.H,obj.n);
        [X,Y]=meshgrid(x,y);
        num=obj.BINT((UP-U).^2);
        denom=obj.BINT((U).^2);
        M=(UP-U).^2;
        numM=obj.BINT(M(nB:end-nB,nB:end-nB));
        denomM=obj.BINT(U(nB:end-nB,nB:end-nB).^2);


    end



        function []=DiagBcomp(obj)
            disp(" Computing DB for "+ obj.Name)
             % generate the kinetic energy operator in finite difference form
             h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
             q_el = 1.60218e-19; % elementary charge, in coulombs
             Mass=9.10938215e-31; % Rest mass of particle in k
             
             %SnapShots=size(obj.Modes,3);
             obj.DB1=zeros(obj.LSize);
             obj.DB2=zeros(obj.LSize);
              for i=1:obj.LSize
                 etai=obj.Modes(:,:,i);
                [etai_dx,etai_dy]=gradient(obj.Modes(:,:,i),obj.dx,obj.dy);
                for j=i:obj.LSize
                    etaj=obj.Modes(:,:,j);
                    [etaj_dx,etaj_dy]=gradient(obj.Modes(:,:,j),obj.dx,obj.dy);
                    
                    %integrate top face
                    FT=obj.gT(etai_dy.*etaj+etai.*etaj_dy);
                    I=obj.xint(FT);

                    %Bottom Face
                    FB=obj.gB(etai_dy.*etaj+etai.*etaj_dy);
                    I=I-obj.xint(FB); % minus accounts for normal

                    % Compute Right Face
                    FR=obj.gR(etai_dx.*etaj+etai.*etaj_dx);
                    I=I+obj.yint(FR);

                    % Compute left face
                    FL=obj.gL(etai_dx.*etaj+etai.*etaj_dx);
                    I=I-obj.yint(FL); % minus accounts for normal

                    I=I*-1/4*h_bar^2/(Mass*q_el);

                    I2=obj.xint(obj.gT(etai.*etaj));
                    I2=I2+obj.xint(obj.gB(etai.*etaj));
                    I2=I2+obj.yint(obj.gR(etai.*etaj));
                    I2=I2+obj.yint(obj.gL(etai.*etaj));
                    I2=I2*(h_bar.^2/(2*Mass*q_el)); % scale
                    obj.DB1(i,j)=I;
                    obj.DB2(i,j)=I2;
                end
              end
              for i=1:obj.LSize
                   obj.DB1(i,i)=obj.DB1(i,i)/2;
                   obj.DB2(i,i)=obj.DB2(i,i)/2;
              end
              obj.DB1=obj.DB1'+obj.DB1;
              obj.DB2=obj.DB2'+obj.DB2;

        end
        function A=ComputeA(obj)
            %first compute A matrix
            SnapShots=size(obj.WFData,3); % get number of snapshots
            A=zeros(SnapShots); % create the A matrix
            for i=1:SnapShots
                for j=i:SnapShots
                    Product=obj.WFData(:,:,i).*obj.WFData(:,:,j);
                    A(i,j)=obj.BINT(Product);
                end
            end
            for i=1:SnapShots
                A(i,i)=A(i,i)/2;
            end
             A=A'+A;
        end
 function []=ComputeModes(obj)
            x=linspace(0,obj.W,obj.n);
            y=linspace(0,obj.H,obj.n);
            obj.x=x;
            obj.y=y;
            obj.dx=x(2)-x(1);
            obj.dy=y(2)-y(1);
            obj.dA=obj.dx*obj.dy;
            SnapShots=size(obj.WFData,3); % get number of snapshots
            disp("Computing A Matrix for "+ obj.Name);
            A=obj.ComputeA();
            [U,Lambda]=eig(A);
            Lambda=diag(Lambda);
            [Lambda,I]=sort(Lambda,"descend");
            obj.EigenValues=Lambda/Lambda(1);
            U=U(:,I); % sort the u vectors
            % plot the eigenvalues
            figure;
            semilogy(1:SnapShots,Lambda/Lambda(1),'linewidth',3,"DisplayName",obj.Name);
            xlabel("Mode Index","FontSize",16)
            ylabel("Eigenvalue","FontSize",16)
            legend()
            [status, msg, msgID]=mkdir("./Library/"+obj.Name);
            savefig("./Library/"+obj.Name) % save the figure to the library
            % now compute the modes
            obj.Modes=zeros(size(obj.WFData,1),size(obj.WFData,2),obj.LSize);
            obj.ModesMult=zeros(obj.n*obj.n,obj.LSize);
            disp("Creating POD Modes for "+obj.Name)
            MatData=reshape(obj.WFData,obj.n^2,SnapShots);
            obj.TM=zeros(1,obj.n,obj.LSize);
            obj.BM=zeros(1,obj.n,obj.LSize);
            obj.LM=zeros(obj.n,1,obj.LSize);
            obj.RM=zeros(obj.n,1,obj.LSize);
            obj.TG=zeros(1,obj.n,obj.LSize);
            obj.BG=zeros(1,obj.n,obj.LSize);
            obj.LG=zeros(obj.n,1,obj.LSize);
            obj.RG=zeros(obj.n,1,obj.LSize);

            for i=1:obj.LSize               
               obj.Modes(:,:,i)=reshape(MatData*U(:,i),obj.n,obj.n); 

                % normalize the POD modes
                obj.Modes(:,:,i)=normalize(obj.Modes(:,:,i),obj.dA).*obj.Modes(:,:,i);
                obj.ModesMult(:,i)=reshape(obj.Modes(:,:,i),obj.n*obj.n,1);
                obj.TM(:,:,i)=obj.gT(obj.Modes(:,:,i));
                obj.BM(:,:,i)=obj.gB(obj.Modes(:,:,i));
                obj.LM(:,:,i)=obj.gL(obj.Modes(:,:,i));
                obj.RM(:,:,i)=obj.gR(obj.Modes(:,:,i));
            end

            obj.FX=zeros(size(obj.WFData,1),size(obj.WFData,2),SnapShots);
            obj.FY=zeros(size(obj.WFData,1),size(obj.WFData,2),SnapShots);
            for i=1:obj.LSize
                [FX,FY]=gradient(obj.Modes(:,:,i),obj.dx,obj.dy);
                obj.FX(:,:,i)=FX;
                obj.FY(:,:,i)=FY;
                obj.TG(:,:,i)=obj.gT(obj.FY(:,:,i));
                obj.BG(:,:,i)=obj.gB(obj.FY(:,:,i));
                obj.LG(:,:,i)=obj.gL(obj.FX(:,:,i));
                obj.RG(:,:,i)=obj.gR(obj.FX(:,:,i));
            end

        end
       
        % function []=ComputeModes(obj)
        %     x=linspace(0,obj.W,obj.n);
        %     y=linspace(0,obj.H,obj.n);
        %     obj.x=x;
        %     obj.y=y;
        %     obj.dx=x(2)-x(1);
        %     obj.dy=y(2)-y(1);
        %     obj.dA=obj.dx*obj.dy;
        %     SnapShots=size(obj.WFData,3); % get number of snapshots
        %     disp("Computing SVD Matrix for "+ obj.Name);
        %     [U,Lambda,~]=svd(reshape(obj.WFData,obj.n*obj.n,SnapShots));
        %     Lambda=diag(Lambda).^2;
        %     % A=obj.ComputeA();
        %     % [U,Lambda]=eig(A);
        %     % Lambda=diag(Lambda);
        %     % [Lambda,I]=sort(Lambda,"descend");
        %     obj.EigenValues=Lambda/Lambda(1);
        % 
        %     % plot the eigenvalues
        %     figure;
        %     semilogy(1:1:SnapShots,Lambda/Lambda(1),'linewidth',3,"DisplayName",obj.Name);
        %     xlabel("Mode Index","FontSize",16)
        %     ylabel("Eigenvalue","FontSize",16)
        %     legend()
        %     [status, msg, msgID]=mkdir("./Library/"+obj.Name);
        %     savefig("./Library/"+obj.Name) % save the figure to the library
        %     % now compute the modes
        %     obj.ModesMult=zeros(obj.n*obj.n,obj.LSize);
        %     obj.Modes=zeros(size(obj.WFData,1),size(obj.WFData,2),obj.LSize);
        %     disp("Creating POD Modes for "+obj.Name)
        % 
        %     obj.TM=zeros(1,obj.n,obj.LSize);
        %     obj.BM=zeros(1,obj.n,obj.LSize);
        %     obj.LM=zeros(obj.n,1,obj.LSize);
        %     obj.RM=zeros(obj.n,1,obj.LSize);
        %     obj.TG=zeros(1,obj.n,obj.LSize);
        %     obj.BG=zeros(1,obj.n,obj.LSize);
        %     obj.LG=zeros(obj.n,1,obj.LSize);
        %     obj.RG=zeros(obj.n,1,obj.LSize);
        % 
        %     for i=1:obj.LSize
        % 
        %        obj.Modes(:,:,i)=reshape(U(:,i),obj.n,obj.n);              
        %         % normalize the POD modes
        %         obj.Modes(:,:,i)=normalize(obj.Modes(:,:,i),obj.dA).*obj.Modes(:,:,i);
        %         obj.ModesMult(:,i)=reshape(obj.Modes(:,:,i),obj.n*obj.n,1);
        %         obj.TM(:,:,i)=obj.gT(obj.Modes(:,:,i));
        %         obj.BM(:,:,i)=obj.gB(obj.Modes(:,:,i));
        %         obj.LM(:,:,i)=obj.gL(obj.Modes(:,:,i));
        %         obj.RM(:,:,i)=obj.gR(obj.Modes(:,:,i));
        %     end
        % 
        %     obj.FX=zeros(size(obj.WFData,1),size(obj.WFData,2),obj.LSize);
        %     obj.FY=zeros(size(obj.WFData,1),size(obj.WFData,2),obj.LSize);
        %     for i=1:obj.LSize
        %         [FX,FY]=gradient(obj.Modes(:,:,i),obj.dx,obj.dy);
        %         obj.FX(:,:,i)=FX;
        %         obj.FY(:,:,i)=FY;
        %         obj.TG(:,:,i)=obj.gT(obj.FY(:,:,i));
        %         obj.BG(:,:,i)=obj.gB(obj.FY(:,:,i));
        %         obj.LG(:,:,i)=obj.gL(obj.FX(:,:,i));
        %         obj.RG(:,:,i)=obj.gR(obj.FX(:,:,i));
        %     end
        % 
        % end

        function A=ComputeAU(obj)
            %first compute A matrix
            SnapShots=size(obj.UData,3); % get number of snapshots
            A=zeros(SnapShots); % create the A matrix
            for i=1:SnapShots
                for j=i:SnapShots
                    Product=obj.UData(:,:,i).*obj.UData(:,:,j);
                    A(i,j)=obj.BINT(Product);
                end
            end
            for i=1:SnapShots
                A(i,i)=A(i,i)/2;
            end
             A=A'+A;
        end

    % function []=ComputeUModes(obj,USize)
    %         obj.NumUM=USize;
    %         SnapShots=size(obj.UData,3); % get number of snapshots
    %         disp("Computing svd Matrix for "+ obj.Name);
    %         [U,Lambda,V]=svd(reshape(obj.UData,obj.n.^2,SnapShots));
    %         Lambda=diag(Lambda).^2
    %         Lambda=Lambda/Lambda(1);
    %         obj.EigenValuesU=Lambda;
    % 
    % 
    % 
    %         % plot the eigenvalues
    %         plt=figure;
    %         semilogy(1:SnapShots,Lambda/Lambda(1),'linewidth',3,"DisplayName",obj.Name);
    %         xlabel("Mode Index","FontSize",16)
    %         ylabel("Eigenvalue for Potential data","FontSize",16)
    %         savefig(plt,"./Library/U_Spectrum"+obj.Name) % save the figure to the library
    %         % now compute the modes
    %         obj.UModes=zeros(size(obj.UData,1),size(obj.UData,2),obj.NumUM);
    %         disp("Creating POD U Modes for "+obj.Name)
    % 
    %         for i=1:obj.NumUM               
    %            obj.UModes(:,:,i)=reshape(U(:,i),obj.n,obj.n);              
    %         end
    % 
    %     end
      function []=ComputeUModes(obj,USize)
            SnapShots=size(obj.UData,3); % get number of snapshots
            disp("Computing A_U Matrix for "+ obj.Name);
            A=obj.ComputeAU();
            [U,Lambda]=eig(A);
            Lambda=diag(Lambda);
            [Lambda,I]=sort(Lambda,"descend");
            Lambda=Lambda/Lambda(1);
            obj.EigenValuesU=Lambda;
            obj.NumUM=USize;

            U=U(:,I); % sort the u vectors
            % plot the eigenvalues
            figure;
            semilogy(1:SnapShots,Lambda/Lambda(1),'linewidth',3,"DisplayName",obj.Name);
            xlabel("Mode Index","FontSize",16)
            ylabel("Eigenvalue for Potential data","FontSize",16)
            savefig("./Library/U_Spectrum"+obj.Name) % save the figure to the library
            % now compute the modes
            obj.UModes=zeros(size(obj.UData,1),size(obj.UData,2),obj.NumUM);
            disp("Creating POD U Modes for "+obj.Name)
            MatData=reshape(obj.UData,obj.n^2,SnapShots);

            for i=1:obj.NumUM               
               obj.UModes(:,:,i)=reshape(MatData*U(:,i),obj.n,obj.n);              
                % normalize the POD modes
               obj.UModes(:,:,i)=normalize(obj.UModes(:,:,i),obj.dA).*obj.UModes(:,:,i);
            end

        end


        function [] =DebugM(obj)
            x=linspace(0,obj.W,obj.n);
            y=linspace(0,obj.H,obj.n);
            [X,Y]=meshgrid(x,y);
            for i=1:10
                figure;
                contourf(X,Y,obj.Modes(:,:,i));
            end
        end


    end

end