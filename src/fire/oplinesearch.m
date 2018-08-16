function [Jmin,smin,um,Jlow] = oplinesearch(u,f,R,dir,p)
% [Jmin,umin,um,step_opt,Jlow] = oplinesearch(u,f,dir,p,step)
%inputs:
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
%outputs:
%   Jmin        minim value of J in the patch
%   smin        step where J is minimum in the patch
%   um          new array doing descent gradient methot from u
%   Jlow        value of J in the patch before line search

[m,n]=size(u);
step_low=0;
step_high=p.max_step;
Jmins=zeros(p.max_depth,1);
Jdiffs=zeros(p.max_depth,1);
smins=zeros(p.max_depth,1);
[Jlow,~]=f(u,R,p);
for d=1:p.max_depth
    % defining the steps and the line search array Jls
    steps=linspace(step_low,step_high,p.nmesh+2);
    Jls=Jlow;
    for i=2:p.nmesh+2
        % new u
        v=u+steps(i)*dir;
        %err
        %v=min(v,p.umax(1));
        %v=max(v,p.umin(1));
        % difference on f
        [fv,~]=f(v,R,p);
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
            K=10;   
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
%um=u+smin*(dir-u);
um=u+smin*dir;
%um=min(um,p.umax(1));
%um=max(um,p.umin(1));
%{
function err
    fprintf('oplinesearch: constraint relative error %d, scal %d\n',norm(H*v(:)-g)/norm(g),norm(g));
end
%}
end

