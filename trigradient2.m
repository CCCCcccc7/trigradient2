function [zx, zy, zxx, zyy, zxy] = trigradient2(X,Y,Z,T,M)
%% Definition
% trigradient2 computes first and second order derivatives of a function
%   on a 2D-points region.   
% 
%% Usage
% Input: 
%   X= Vector with x-coordinates.
%   Y= Vector with y-coordinates.
%   Z= Matrix with function values on each point. If Z has multiple
%       columns the derivatives for each column are calculated.
%  Optional argument:
%   T= Triangulation (Nx3 matrix with polygon vertices). If not given the delaunay triangulation of X,Y is used.
%   M= Method used for calculation. Default value is 0.
%      0: Fully vectorized approach with one large equation system. Fast. 
%      1: Multiple small equations systems. Slower but depending on the input values more accurate.
%
% Output:
%   zx=dz/dx
%   zy=dz/dy
%   zxx=d^2z/dx^2
%   zxy=d^2z/(dx dy)
%   zyy=d^2z/dy^2
%
%% Method
% The derivatives are calculated with a least squares linear regression. The system of
%   equations is set up with taylor series from each point to the adjacent 
%   vertices. If a vertex is connected to less than 5 vertices then vertices in 
%   2-edge distance are used too.
%
%% Example:
% [x y]=meshgrid(-1:.1:1,-1:.1:1);
% x=x(:);
% y=y(:);
% n=1:length(x);
% dt=delaunay(x,y);
% f=exp((x.^2+y.^2)/100).*sin(x).*cos(y);
% dfdx=exp((x.^2+y.^2)/100).*cos(y).*(cos(x)+sin(x).*x/50);
% dfdy=exp((x.^2+y.^2)/100).*sin(x).*(-sin(y)+cos(y).*y/50);
% d2fdx2=exp((x.^2+y.^2)/100).*cos(y).*(sin(x).*x.^2/2500+cos(x).*x/25-sin(x).*49/50);
% d2fdy2=exp((x.^2+y.^2)/100).*sin(x).*(cos(y).*y.^2/2500-sin(y).*y/25-cos(y).*49/50);
% d2fdxy=exp((x.^2+y.^2)/100).*((y.*cos(y)/50).*(cos(x)+x.*sin(x)/50)-sin(y).*(cos(x)+x.*sin(x)/50));
% [dfdxE, dfdyE, d2fdx2E, d2fdy2E, d2fdxyE]=trigradient2(x,y,f,dt);
% figure
%  subplot(5,1,1)
%  plot(n,dfdx,n,dfdxE)
%  ylabel('df/dx')
%  legend('Exact value','trigradient2.m approx.','Location','EastOutside')
%  subplot(5,1,2)
%  plot(n,dfdy,n,dfdyE)
%  ylabel('df/dy')
%  subplot(5,1,3)
%  plot(n,d2fdx2,n,d2fdx2E)
%  ylabel('d^2f/dx^2')
%  subplot(5,1,4)
%  plot(n,d2fdy2,n,d2fdy2E)
%  ylabel('d^2f/dy^2')
%  subplot(5,1,5)
%  plot(n,d2fdxy,n,d2fdxyE)
%  ylabel('d^2f/(dx dy)')
%  xlabel('Point number')
%
%% See also
% trigradient by Mick Warehime.
%
%% Author
% Gypaets
%
%% Function
% Default values
if nargin==4
    M=0;
elseif nargin==3
    T=delaunay(X,Y);
    M=0;
end

v=version;
v=str2num(v(end-5:end-2)); % Matlab version year

nV=length(X); % Number of vertices
nT=length(T); % Number of polygons
nZ=length(Z(1,:)); % Number of columns in function
nvpK=5; % Number of variables to be calculated per point

T=sort(T,2);
XY=[X Y];

% Adjacent vertices
[lr lc]=find(sparse(reshape(T(:,[1 1 2]).',[],1),reshape(T(:,[2 3 3]).',[],1),1,nV,nV));
adjacencyMatrix=sparse(lr,lc,1,nV,nV)+sparse(lr,lc,1,nV,nV).';

% Vertices connected to 4 or less points
sparseV=find(sum(adjacencyMatrix)'<=4);

% 2-polygon distance adjacency matrix
[lr lc]=find(adjacencyMatrix+adjacencyMatrix^2);
adjacencyMatrix2=sparse(lr,lc,1,nV,nV)-sparse(1:nV,1:nV,1);

% For derivative calculation used adjacency matrix
adjacencyMatrix(:,sparseV)=adjacencyMatrix2(:,sparseV);

% Set up and solve equation system
if M==0
% Fully vectorized method: one large equation system
    % LHS matrix
    rowsNumber=sum(adjacencyMatrix(:));
    
    if v >= 2015
        columnInd=repelem(1:nvpK:nvpK*nV,sum(adjacencyMatrix)');
    elseif v < 2015
        clens = cumsum(sum(adjacencyMatrix));
        idx(clens(end))=0;
        idx([1 clens(1:end-1)+1]) = diff([0 1:nvpK:nvpK*nV]);
        columnInd = cumsum(idx);
    end
    allColumnInd=[columnInd, columnInd+1, columnInd+2, columnInd+3, columnInd+4];
    
    [lr lc]=find(adjacencyMatrix);
    dX=X(lr)-X(lc);
    dY=Y(lr)-Y(lc);
    
    lhsMatrix=sparse(repmat(1:rowsNumber,1,nvpK),allColumnInd,[dX; dY; 0.5*dX.^2; 0.5*dY.^2; dX.*dY],rowsNumber,nvpK*nV);
    
    % RHS 'vector'
    rhsVector=Z(lr,:)-Z(lc,:);
    
    % Derivative calculation
    dZ=lhsMatrix\rhsVector;
    
    % Assign values
    zx=dZ(1:5:end,:);
    zy=dZ(2:5:end,:);
    zxx=dZ(3:5:end,:);
    zyy=dZ(4:5:end,:);
    zxy=dZ(5:5:end,:);

elseif M==1
% Multiple small equation systems
    % Preallocate result matrix
    dZ{nV}=zeros([nZ 5]);
    
    % Loop over all mesh vertices
    for lk=1:nV
        laK=find(adjacencyMatrix(:,lk)); % Vertices used for calculation on vertex lk.
        
        % Loop over neighbouring vertices to compute LHS Matrix and RHS vector
        clear lA lbV
        for li=1:length(laK)
            lxy=XY(laK(li),:)-XY(lk,:); % Distance vector to neighbour vertex
            lA(li,:)=[lxy lxy.^2/2 prod(lxy)]; % LHS Matrix column
            lbV(li,:)=Z(laK(li),:)-Z(lk,:); % RHS Vector column
        end
        
        % Derivative calculation
        if rank(lA)>=5
            % Least square fitting
            dZ{lk}=(lA\lbV)';
        elseif rank(lA)<=4
            % Display warning and set derivatives to zero if LHS Matrix has not full rank.
            dZ{lk}=zeros([nZ,5]);
            disp(['trigradient2.m: regression matrix on vertex ' num2str(lk) ' has a rank < 5. Derivatives at this vertex have been set to 0.']);
        end
    end

    % Assign values
    dZm=permute(cat(3,dZ{:}),[3 1 2]);
    zx=dZm(:,:,1); 
    zy=dZm(:,:,2);
    zxx=dZm(:,:,3);
    zyy=dZm(:,:,4);
    zxy=dZm(:,:,5);
 end
