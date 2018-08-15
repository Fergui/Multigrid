function H = interop_bary(X,Y,xq,yq)
% Call:
% function H = interop_bary(X,Y,xq,yq)
%
% Description:
% Compute the interpolation operator from rectangular grid of 
% points with coordinates X and Y to points with coordinates xq and yq
% by the method of barycentric coordinates
%
% Inputs:
%      X,Y        matrices of grid cooredinates
%      xq,yq      coordinates of points to interpolate to
% Outputs:
%      H          Sparse matrix, interpolation operator from values at
%                 X, Y ordered by columns first
%
%-------------------------------------------------------------------------

[m,n]=size(X);  % grid size
nq=numel(xq);  % number of points to interpolate to

% indices of triangles in the grid
%
%      (i,j)----(i,j+1)
%      |       / |
%      | 1   /   |
%      |   /  2  |
%    (i+1,j)---- (i+1,j+1)

% to compute barycentric coordinates b1 b2 b3 of point x,y in triangle
% with vertices [xv1,yv1], [xv2,yv2], [xv3,yv3]:
%  xv1*b1 + xv2*b2 + xv3*b3 = x
%  yv1*b1 + yv2*b2 + yv3*b3 = y
%      b1 +     b2 +     b3 = 1

testing = 0;

[xs,xi]=sort(xq);  % sort the coordinates
[ys,yi]=sort(yq);
H = sparse([],[],[],nq,m*n,nq*3);
if testing,
    H_orig = sparse([],[],[],nq,m*n,nq*3);
end
A=ones(3,3);
C=[xq(:)';yq(:)';ones(1,nq)];  % right hand side 
[sxq,ixq]=sort(xq); % for computing shapefile points within bounds
[syq,iyq]=sort(yq);
zq=zeros(nq,1);  % for calculating intersection of indices
for it=1:2
    for i=1:m-1
        for j=1:n-1;
            % indices of triangle vertices in column first numbering
            % node (i,j) stored at location i+m*(j-1) 
            if it==1,
                %  (i,j), (i+1,j), (i,j+1)
                ixv=[i+(j-1)*m,i+1+(j-1)*m,i+j*m];
            else
                %  (i+1,j), (i+1,j+1), (i,j+1)
                ixv=[i+1+(j-1)*m,i+1+j*m,i+j*m];
            end
            % the matrix of the system for barycentric coordinates
            xv=X(ixv);
            yv=Y(ixv);
            A(1,:)=xv;
            A(2,:)=yv;
            % A(3,:)=ones(1,3); % already is
            six=binsearch_bounds(sxq,min(xv),max(xv)); % indices in sorted between
            siy=binsearch_bounds(syq,min(yv),max(yv));
            jjxx=[];
            if ~(isempty(six) | isempty(siy))
                ix=ixq(six);  % selected indices in original numbering
                iy=iyq(siy);
                zq(ix)=ix;     % start on intersection
                in=zq(iy);  % intersection of indices and some zeros
                zq(ix)=0;      % cleanup
                intersection=in(find(in)); % get rid of zeros
                if ~isempty(intersection),
                    CC=C(:,intersection);
                    BB = A\CC;  % right hand side remains the same
                    % B(ii,jj) is interpolation coefficient from node ixv(ii) to point jj
                    % if point jj is in the triangle, which is when all B(:,jj) >=0
                    jjxx = find(~any(BB<0,1)); % testing for <0 is faster than for >=
                    intx = intersection(jjxx);
                end
            end
            if ~isempty(jjxx)  % the triangle contains any points to interpolate to
               H(intx,ixv)=BB(:,jjxx)';
            end        
            if testing,
                % original code
                % the matrix of the system for barycentric coordinates
                A(1,:)=X(ixv);
                A(2,:)=Y(ixv);
                % A(3,:)=ones(1,3); % already is
                B = A\C;  % right hand side remains the same
                % B(ii,jj) is interpolation coefficient from node ixv(ii) to point jj
                % if point jj is in the triangle, which is when all B(:,jj) >=0
                jjx = find(~any(B<0,1)); % testing for <0 is faster than for >=
                if ~isempty(jjx)  % the triangle contains any points to interpolate to
                   H_orig(jjx,ixv)=B(:,jjx)';
                end
            end
        end
    end
end
if testing,
    err=norm(abs(H_orig-H),1)
end
end
