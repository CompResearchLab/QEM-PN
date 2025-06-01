classdef SingleGauss <Block & handle % inherite the block and handle class
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties(Access=private)
        SS % shift scale (The maximum shift that is allowed in each direction)
        RMax % max radius
        Rm % min radiuys
        Depth
        
    end

    methods
        function obj = SingleGauss(H,W,SS,RM,Rm,Depth,n)
            %Set height and width, max shift and radius max and min
            obj.H=H; % set height
            obj.W=W; % set width
            obj.SS=SS; % set shift scale 
            obj.RMax=RM; % set radius max
            obj.Rm=Rm; % set radius min
            obj.Depth=Depth; % specify the depth of the gaussian potential
            obj.n=n;
        end
       
        function U= genUR(obj,X,Y,x0,y0)
                U=zeros(size(X,1),size(Y,2));
                EW=X(end); % get entire width
                EL=Y(end); % get entire legnth
                XS=[0, 0  ,   0 , EW , EW, EW ,-EW , -EW , -EW];
                YS=[0, EL , -EL , EL ,  0, -EL, EL ,   0 ,  -EL];
                sx=(obj.SS+obj.SS)*rand()-obj.SS; % generate a random shift x
                sy=(obj.SS+obj.SS)*rand()-obj.SS; % generate a random shift y
                Radius=(obj.RMax-obj.Rm)*rand()+obj.Rm;
                for i=1:length(XS)
                    cX=obj.W/2+x0; % get the x center of the block
                    cY=obj.H/2+y0; % get the y center of the block
                    XN=(X-cX-sx-XS(i)); % shifted cordinates X
                    YN=(Y-cY-sy-YS(i)); % shifted cordinates Y
                    U=U-obj.Depth*exp(-XN.^2/(2*Radius.^2)-YN.^2/(2*Radius.^2));
                end
        end


        % vars argument is described as follows:
        % first parameter is shift x 
        % second parameter is shift y
        % Third parameter is the Radius of the sinc function
        function U=genUT(obj,X,Y,x0,y0,vars)
             U=zeros(size(X,1),size(Y,2));
                EW=X(end); % get entire width
                EL=Y(end); % get entire legnth
                XS=[0, 0  ,   0 , EW , EW, EW ,-EW , -EW , -EW];
                YS=[0, EL , -EL , EL ,  0, -EL, EL ,   0 ,  -EL];
                sx=vars(1); % generate a random shift x
                sy=vars(2); % generate a random shift y
                Radius=vars(3);
                for i=1:length(XS)
                    cX=obj.W/2+x0; % get the x center of the block
                    cY=obj.H/2+y0; % get the y center of the block
                    XN=(X-cX-sx-XS(i)); % shifted cordinates X
                    YN=(Y-cY-sy-YS(i)); % shifted cordinates Y
                    U=U-obj.Depth*exp(-XN.^2/(2*Radius.^2)-YN.^2/(2*Radius.^2));
                end
           
        
        end
    end
end