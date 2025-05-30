function B=FDTP(X,Y,n)
        % generate the kinetic energy operator in finite difference form
        h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
        q_el = 1.60218e-19; % elementary charge, in coulombs
        Mass=9.10938215e-31; % Rest mass of particle in kg
        %------------------------Convert Input spacial variables meters-----------%
        % find x and y cordinates
        x=X(1,:);
        dh=x(2)-x(1);
         % Form the eigenvalue problem
        Index=zeros(3,(n-1)*9); % create a sparse matrix
        %G=sparse((n-1)*(n-1));
        NeiC1=-h_bar^2/(2*Mass*q_el)*16/(12*dh^2); % get neighbor coefficents
        NeiC2=-h_bar^2/(2*Mass*q_el)*(-1)/(12*dh^2);
        CP=-h_bar^2/(2*Mass*q_el)*(-60)/(12*dh^2);
        I=1;
        for j=1:n-1 %column index
            for i=1:n-1 %row index
                GI=GlobalI(n-1,n-1,i,j);
                CNeiC1=NeiC1;
                CNeiC2=NeiC2;
                CCP=CP;
                Index(:,I)=[GI;GI;CCP]; % put in center value
                I=I+1;
                % put in horizontal neighbors
                if j==1 % left boundary
                    % put in left neighbors
                    NI1=GlobalI(n-1,n-1,i,n-1);% get index of first left boundary
                    NI2=GlobalI(n-1,n-1,i,n-2); % get index of second left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in right neighbors
                    NI1=GlobalI(n-1,n-1,i,j+1); % get index of first right neighbor
                    NI2=GlobalI(n-1,n-1,i,j+2); % get index of second right neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                elseif j==2
                    % put in left neighbors
                    NI1=GlobalI(n-1,n-1,i,1);% get index of first left boundary
                    NI2=GlobalI(n-1,n-1,i,n-1); % get index of second left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in right neighbors
                    NI1=GlobalI(n-1,n-1,i,j+1); % get index of first right neighbor
                    NI2=GlobalI(n-1,n-1,i,j+2); % get index of second right neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                elseif j==(n-1)% Right boundary
                    % put in right neibors
                    NI1=GlobalI(n-1,n-1,i,1); % get index for first right neighbor
                    NI2=GlobalI(n-1,n-1,i,2); % get index for second right neighbor
                    Index(:,I)=[GI,NI1,CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI,NI2,CNeiC2];
                    I=I+1;
                    % put in left neighbors
                    NI1=GlobalI(n-1,n-1,i,j-1);
                    NI2=GlobalI(n-1,n-1,i,j-2);
                    Index(:,I)=[GI,NI1,CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI,NI2,CNeiC2];
                    I=I+1;
                elseif j==(n-2)
                     % put in right neibors
                    NI1=GlobalI(n-1,n-1,i,n-1); % get index for first right neighbor
                    NI2=GlobalI(n-1,n-1,i,1); % get index for second right neighbor
                    Index(:,I)=[GI,NI1,CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI,NI2,CNeiC2];
                    I=I+1;
                    % put in left neighbors
                    NI1=GlobalI(n-1,n-1,i,j-1);
                    NI2=GlobalI(n-1,n-1,i,j-2);
                    Index(:,I)=[GI,NI1,CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI,NI2,CNeiC2];
                    I=I+1;
                else
                     % put in left neighbors
                    NI1=GlobalI(n-1,n-1,i,j-1);% get index of first left boundary
                    NI2=GlobalI(n-1,n-1,i,j-2); % get index of second left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in right neighbors
                    NI1=GlobalI(n-1,n-1,i,j+1); % get index of first right neighbor
                    NI2=GlobalI(n-1,n-1,i,j+2); % get index of second right neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                end
                % put in vertical neighbors
                if i==1
                    % put in bottom neighbors
                    NI1=GlobalI(n-1,n-1,n-1,j);% get index of first bottom boundary
                    NI2=GlobalI(n-1,n-1,n-2,j); % get index of bottom left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in top neighbors
                    NI1=GlobalI(n-1,n-1,i+1,j); % get index of first top neighbor
                    NI2=GlobalI(n-1,n-1,i+2,j); % get index of second top neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                elseif i==2
                     % put in bottom neighbors
                    NI1=GlobalI(n-1,n-1,1,j);% get index of first bottom boundary
                    NI2=GlobalI(n-1,n-1,n-1,j); % get index of bottom left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in top neighbors
                    NI1=GlobalI(n-1,n-1,3,j); % get index of first top neighbor
                    NI2=GlobalI(n-1,n-1,4,j); % get index of second top neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                elseif i==n-1
                     % put in bottom neighbors
                    NI1=GlobalI(n-1,n-1,i-1,j);% get index of first bottom boundary
                    NI2=GlobalI(n-1,n-1,i-2,j); % get index of bottom left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in top neighbors
                    NI1=GlobalI(n-1,n-1,1,j); % get index of first top neighbor
                    NI2=GlobalI(n-1,n-1,2,j); % get index of second top neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                elseif i==n-2
                 % put in bottom neighbors
                    NI1=GlobalI(n-1,n-1,i-1,j);% get index of first bottom boundary
                    NI2=GlobalI(n-1,n-1,i-2,j); % get index of bottom left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in top neighbors
                    NI1=GlobalI(n-1,n-1,n-1,j); % get index of first top neighbor
                    NI2=GlobalI(n-1,n-1,1,j); % get index of second top neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                else
                    % put in bottom neighbors
                    NI1=GlobalI(n-1,n-1,i-1,j);% get index of first bottom boundary
                    NI2=GlobalI(n-1,n-1,i-2,j); % get index of bottom left boundary
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2]; 
                    I=I+1;
                    % put in top neighbors
                    NI1=GlobalI(n-1,n-1,i+1,j); % get index of first top neighbor
                    NI2=GlobalI(n-1,n-1,i+2,j); % get index of second top neighbor
                    Index(:,I)=[GI;NI1;CNeiC1];
                    I=I+1;
                    Index(:,I)=[GI;NI2;CNeiC2];
                    I=I+1;
                end
            end
        end
     B=sparse(Index(1,:),Index(2,:),Index(3,:),(n-1)*(n-1),(n-1)*(n-1));
end

function [I]=GlobalI(TR,RC,CR,CC)
    I=(CC-1)*TR+CR;
end