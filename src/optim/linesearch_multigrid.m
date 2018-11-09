function [Jmin,smin,um,Jlow,Cm,pm] = linesearch_multigrid(u,f,R,dir,vmask,p)
% Call:
% [Jmin,smin,um,Jlow,Cm] = linesearch_multigrid(u,f,R,dir,p)
%
% Description:
% Exact line search in the specify direction
%
% Inputs:
%   u        starting array
%   f        function to make the line search
%   dir      direction to search
%   vmask    mask where to optimize u
%   p        structure with:
%      nmesh       number of mesh points in each search
%      max_depth   max number of searchs
%      min_depth   min number of searchs
%      umin        array, minimal value of u
%      umax        array, maximal value of u
%      H           matrix, interpolation operator in the patch
%      g           array, right side vector of constraints Hu=g
% Outputs:
%   Jmin        minim value of J in the patch
%   smin        step where J is minimum in the patch
%   um          new array doing descent gradient methot from u
%   Jlow        value of J in the patch before line search
%   Cm          contribution matrix where J is minimum
%   pm          penalty value where J is minimum
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
% Modify from another version by Jan Mandel
%-------------------------------------------------------------------------

%figure(10), subplot(2,2,1), mesh(u), view([0,1]), title('First guess of line search'), drawnow;
step_low=0;
step_high=p.max_step;
Jmins=zeros(p.max_depth,1);
Jdiffs=zeros(p.max_depth,1);
smins=zeros(p.max_depth,1);
pmin=zeros(p.max_depth,1);
Cs=cell(p.max_depth,1);
p.vmask=vmask;
[C0,Jlow]=f(u,R,p);
[cm,cn]=size(C0);
for d=1:p.max_depth
    % defining the steps and the line search array Jls
    steps=linspace(step_low,step_high,p.nmesh+2);
    Jls=Jlow;
    ps=zeros(size(steps));
    Cn=C0(:);
    for i=2:p.nmesh+2
        % new u
        v=u+steps(i)*dir;
        %subplot(2,2,2), mesh(v), view([0,1]), title('New solution'), drawnow;
        % difference on f
        [C,fv]=f(v,R,p);
        Cn(:,i)=C(:);
        % penalty term
        if p.penalty
            penalty=zeros(size(u));
            penalty(2:end-1,2:end-1)=max(min(v(1:end-2,2:end-1),min(v(3:end,2:end-1),min(v(2:end-1,1:end-2),v(2:end-1,3:end))))-v(2:end-1,2:end-1),0);
            pen=penalty; 
            pen(~p.vmask)=0;
            %scal=min(v(:))-max(v(:));
            %pen=(penalty/(scal+realmin)).^2*abs(Jls(i-1)-fv)/2;
            %subplot(2,2,3), pcolor(penalty), colorbar, drawnow;
            if isfield(p,'K')
                K=p.K;
            else
                K=10;
            end
            % penalty value
            ps(i)=K*sum(pen(:));
            % new function value
            Jls(i)=fv+ps(i);
        else
            % new function value
            Jls(i)=fv;
        end
        
        % conditions of the line search
        out=0;
        step_opt = steps(i);
        if Jls(i)>Jls(i-1)
            step_opt = steps(i-1);
            out=0;
            break
        end
        if (Jls(i)<Jls(1) && d>=p.min_depth)
            out=1;
            break
        end
    end
    % minimizing through all the step sizes
    [Jmins(d),ndx]=min(Jls);
    Jdiffs(d)=Jls(ndx)-Jlow;
    smins(d)=steps(ndx);
    pmin(d)=ps(ndx);
    Cs{d}=reshape(Cn(:,ndx),cm,cn);
    % redifine bounds
    low=max(ndx-1,1);
    high=min(ndx+1,p.nmesh+2);
    Jlow=Jls(low);
    step_low=steps(low);
    if out
        break
    end
    if high<p.nmesh+2
        step_high=steps(high);
    else
        step_high=steps(high)*2;
    end
end
% Taking the best deep level and the best alpha candidate
[~,ndx]=min(Jdiffs);
Jmin=Jmins(ndx);
smin=smins(ndx);
Cm=Cs{ndx};
pm=pmin(ndx);
um=u+smin*dir;
%subplot(2,2,4), mesh(um), view([0,1]), title('Final result'), drawnow;
end