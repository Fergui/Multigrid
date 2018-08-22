function [um,Jglobal]=multigrid(u,p)
% Call:
% um=multigrid(u,p)
%
% Description:
% Computes a linesearch technique using multigrid approach
%
% Inputs:
%   u    Initial approximation
%   p    Matlab structure containing:
%               R           matrix, rate of spread, on same nodes as u
%               H           sparse matrix, interpolation operator for all
%                           the constraints Hu=g
%               g           vector, right side of the constraints (Hu=g)
%               ofunc       matlab function, objective function comparing 
%                           x=||grad u||^2 and y=R^2 such that xy=1
%               q           q norm of the computation of J
%               select      handle to upwinding function
%               max_step    max size of the steps to search
%               nmesh       number of mesh points in each search
%               max_depth   max number of searchs
%               min_depth   min number of searchs
%               umax        array, maximal value of u
%               umin        array, minimal value of u    
%               mcycle      number of complete cycles in the multigrid
%               lm          number of coarse mesh step cycles
%               seq         array, sequence of coarse mesh steps
%               multigrid   array, number of cycles in each mesh step size
%               ros         boolean variable to know if compute or not a dynamic ROS
%               rec         boolean variable to know if record the optimization
%               vmask       (optional) matrix, true where the values of level set function can change
% Outputs:
%   um        Final solution of the multigrid method
%   Jglobal   Array of objective function values after each iteartion
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

% Initialization of the multigrid linesearch with strategy vector
[m,n]=size(u);
lm=size(p.multigrid,2)-1;
seq=lm:-1:0;
um=u;
Jop=[];
Jglobal=[];
r=gJ(um,p);
Jop=[Jop;r];
Jglobal=[Jglobal;r];
% General cycle of multigrid method
for cycle=1:p.mcycle
    % Each cycle will have all the spacings
    for l=1:lm+1
        % Computing the spacing
        s=2^seq(l);
        % Creating the coarse basis function at the given spacing s, phi
        x=(-s:s);
        y=(-s:s);
        [X,Y]=meshgrid(x,y);
        D=(1-abs(X)/s).*(1-abs(Y)/s);
        phi=max(0,D);
        % Plotting the coarse basis function
        subplot(2,2,2)
        if p.ros
            ros=p.R;
            ros(~p.vmask)=nan;
            mesh(ros'), view([0 1]), title('ROS');
        else
            mesh(phi), tit=title(['Bilinear coarse grid function at mesh step ',num2str(s)]);  set(tit,'FontSize',20,'FontWeight','Bold'), axi=zlabel('Fire arrival time'); set(axi,'FontSize',20,'FontWeight','Bold')
        end
        drawnow
        % Generating the possible range nodes (i,j) where we could apply the direction phy
        if isfield(p,'vmask')
            [ii,jj,~]=find(p.vmask);
            mm=min(ii)+s:max(ii)-s;
            nn=min(jj)+s:max(jj)-s;
        else
            mm=1+s:m-s;
            nn=1+s:n-s;
        end
        % Computing the number of iterations on each spacing
        ms=p.multigrid(l);
        % Cycles in each coarse mesh step
        for cs=1:ms
            tic
            fprintf('>>>>> cycle: %d/%d coarse mesh step: %d cycle in coarse mesh step: %d/%d <<<<< \n',cycle,p.mcycle,s,cs,ms);
            for i=mm
                for j=nn
                    if isfield(p,'vmask')
                        if p.vmask(i,j)
                            node=[i,j,s];
                            [um,Jop]=locallinesearch(um,node,phi,Jop,p);
                        end
                    else
                        node=[i,j,s];
                        [um,Jop]=locallinesearch(um,node,phi,Jop,p);
                    end
                end
            end
            r=gJ(um,p);
            Jglobal=[Jglobal;r];
            % Plotting results of each spacing
            subplot(2,2,3)
            plot(0:size(Jglobal,1)-1,Jglobal(1:end),'.-'), tit=title('Objective function after each multigrid iteration'); set(tit,'FontSize',18,'FontWeight','Bold'), axi1=xlabel('Multigrid iteration'); set(axi1,'FontSize',20,'FontWeight','Bold'), axi2=zlabel('Fire arrival time'); set(axi2,'FontSize',20,'FontWeight','Bold'), axi3=ylabel('Objective function value'); set(axi3,'FontSize',20,'FontWeight','Bold')
            subplot(2,2,4)
            ur=um;
            ur(~p.vmask)=nan;
            plot_sol_mesh(p.X,p.Y,ur,p.H,p.g), view([0 1]), tit=title(['T: Cycle=',num2str(cycle),'/',num2str(p.mcycle),' - Coarse mesh step=',num2str(s),' - Cycle in mesh step=',num2str(cs),'/',num2str(ms),' J(T)=',num2str(Jglobal(end))]); set(tit,'FontSize',15,'FontWeight','Bold'), axi=zlabel('Fire arrival time'); set(axi,'FontSize',20,'FontWeight','Bold')
            drawnow
            if p.rec
                record(['multi_',p.exp,'.gif'],p.fig);
            end
            % Updating ROS dynamically
            if p.ros
                p.ds=dyninterp(um,p);
                p.R=ros_file(um,p.ignS,p.ds,p);
                p.R(~p.vmask)=0;
            end
            toc
        end
    end
end
end