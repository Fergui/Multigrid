% Initializing the test
fprintf('Initializing the test...\n');
s=kml2struct('data/doc.kml');
m=500;n=500;
ulon=-106.0347;ulat=36.1729;llon=-106.7677;llat=35.5790;
x=linspace(llon,ulon,m);
y=linspace(llat,ulat,n);
[X,Y]=meshgrid(x,y);
xq=s(1).Lon{3};
xq=xq(1:end-2);
yq=s(1).Lat{3};
yq=yq(1:end-2);

% Add points to the shape
fprintf('Adding points to the shape...\n');
tic
[xqn,yqn]=addshape(X,Y,xq,yq);
toc

% Computing interpolation operator
fprintf('Computing the interpolation operator matrix...\n');
tic
A=interop_bary(X,Y,xqn,yqn);
toc

% Condencing rows of A
fprintf('Condencing rows of the interpolation operator...\n');
tic
[H,rows]=condense(A);
toc

% Tests
B=A(rows,:);
% The same coordinates? Values are different because they are a combination of
% all the points at the same triangle.
B(10,:), H(10,:)
fprintf('# rows of initial H: %d\n',size(A,1));
fprintf('# rows of new H: %d\n',size(H,1));
fprintf('Rank of new H: %d\n',rank(full(H*H')));