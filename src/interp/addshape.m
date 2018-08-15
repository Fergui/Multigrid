function [xqn,yqn] = addshape(X,Y,xq,yq)
%[xqn,yqn] = untitled(X,Y,xq,yq)
% Add shape points if it is necessary depending on the mesh
% in
%      X,Y        matrices of grid cooredinates
%      xq,yq      coordinates of points in the shape
% out
%      xqn,yqn    coordinates of new points in the shape

% Some necessary precomputations 
ns=length(xq);
dx=X(1,2)-X(1,1);
dy=Y(1,2)-Y(1,1);

% Adding points to the shapefile
d=sqrt(([xq(2:end);xq(1)]-xq).^2+([yq(2:end);yq(1)]-yq).^2);
cond=floor(d./sqrt(dx^2+dy^2));
%cond=floor((2*d)./sqrt(dx^2+dy^2));
nc=sum(cond);
xqn=zeros(ns+nc,1); yqn=zeros(ns+nc,1); 
in=0;
for i = 1:ns
    in=in+1;
    xqn(in)=xq(i); 
    yqn(in)=yq(i);
    nci=cond(i);
    for l = 1:nci
        in=in+1;
        if i==ns
            xqn(in)=xq(i)+l*(xq(1)-xq(i))/(nci+1);
            yqn(in)=yq(i)+l*(yq(1)-yq(i))/(nci+1);
        else
            xqn(in)=xq(i)+l*(xq(i+1)-xq(i))/(nci+1); 
            yqn(in)=yq(i)+l*(yq(i+1)-yq(i))/(nci+1);
        end
    end
end
fprintf('The number of points added are %d\n',nc);
%figure, mesh(X,Y,zeros(size(X))), view(2), hold on, plot(xqn,yqn,'+'), hold on, plot(xq,yq,'x'), drawnow


end

