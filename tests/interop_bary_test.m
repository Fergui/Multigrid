function interop_bary_test
% test interop_bary
x=linspace(1,5,5);
y=linspace(1,4,4);
[X,Y]=meshgrid(x,y);
dotest(X,Y,[2:0.25:4],[5:-0.5:1])
dotest(X+0.1*rand(size(X)),Y+0.1*rand(size(Y)),[2:0.25:4],[5:-0.5:1])
end

function dotest(X,Y,xq,yq)
fprintf('interpolation on mesh size %i %i to %i points\n',size(X),numel(xq)); 
H = interop_bary(X,Y,xq,yq);
clf
for i=1:size(H,1)
    fprintf('interpolation to point %g %g :\n',xq(i),yq(i))
    plot(xq(i),yq(i),'*')
    hold on
    jj=find(H(i,:));
    if ~isempty(jj),
        for j=jj
            fprintf('from node %i at %g %g weight %g\n',j,X(j),Y(j),full(H(i,j)))
            % draw the point and the triangle we interpolate from
        end
        fprintf('total weight is %g\n',full(sum(H(i,:))));       
        plot([X(jj),X(jj(1))],[Y(jj),Y(jj(1))],'-')
        plot([X(jj),X(jj(1))],[Y(jj),Y(jj(1))],'o')
    else
        fprintf('the point is not in the mesh envelope\n')
    end
    drawnow
    pause(0.5)
end
press_enter
end

function press_enter
    fprintf('Press enter to continue> \n'),pause
end