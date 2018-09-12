function [Jmin,smin,um,Jlow] = linesearch(u,f,R,dir,p)
% Call:
% [Jmin,smin,um,Jlow] = oplinesearch(u,f,dir,p)
%
% Description:
% Exact line search in the specify direction
%
% Inputs:
%   u        starting array
%   f        function to make the line search
%   dir      direction to search
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
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
% Modify from another version by Jan Mandel
%-------------------------------------------------------------------------

[m,n]=size(u);
step_low=0;
step_high=p.max_step;
Jmins=zeros(p.max_depth,1);
Jdiffs=zeros(p.max_depth,1);
smins=zeros(p.max_depth,1);
Jlow=f(u,R,p);
for d=1:p.max_depth
    % defining the steps and the line search array Jls
    steps=linspace(step_low,step_high,p.nmesh+2);
    Jls=Jlow;
    for i=2:p.nmesh+2
        % new u
        v=u+steps(i)*dir;
        % difference on f
        fv=f(v,R,p);
        % penalty term
        if p.penalty
            penalty=zeros(size(u));
            for k=2:m-1
                for l=2:n-1
                    penalty(k,l)=max(min([v(k-1,l) v(k+1,l) v(k,l-1) v(k,l+1)])-v(k,l),0);
                end
            end
            scal=min(v(:))-max(v(:));
            pen=(penalty/(scal+realmin)).^2*abs(Jls(i-1)-fv)/2;
            K=1000;   
            % new function value
            Jls(i)=fv+K*sum(pen(:));
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
um=u+smin*dir;
end

