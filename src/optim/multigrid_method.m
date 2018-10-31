function mg=multigrid_method(u,p)
tic
[mp,np]=size(u);
sm=findmgfit(mp);
sn=findmgfit(np);
uu=u(sm,sn);
[SM,SN]=meshgrid(sm,sn);
ii=SM'+(SN'-1)*mp;
H=p.H(:,ii(:));
% setup
ratio=2;  
[mm,nn]=size(uu);
% determine number of levels
for levels=1:100,
    mm=(mm-1)/ratio+1; 
    nn=(nn-1)/ratio+1;
    if mm~=round(mm) || nn~=round(nn) || mm<4 || nn<4
        levels
        break
    end
end
[m,n]=size(uu);
% initialize level 1
mg(1).u=uu;
mg(1).m=m;
mg(1).n=n;
mg(1).dx=p.dx;
mg(1).dy=p.dy;
mg(1).nnods=mg(1).m*mg(1).n;
mg(1).H=H;
mg(1).g=p.g;



for i=1:levels
    [p(i).Xc,p(i).Yc]=meshgrid(1:ratio:p(i).m,1:ratio:p(i).n); % coarse indices
    [p(i).C,p(i).I,p(i).Ir] = adini_assemble(p(i).m,p(i).n,adini_local(p(i).dx,p(i).dy)); % problem matrix
    if i==1, p(i).I=[]; end  % prolongation not needed on finest level
    p(i).idx=order_index(p(i).m,p(i).n,2,2); % order of blocks for relaxation
    [p(i).X,p(i).Y]=meshgrid(1:p(i).m,1:p(i).n); % nde indices
    p(i).d=cell(p(i).nnods,1);
    for node=1:p(i).nnods
       ix=3*node+[-2:0];
       p(i).d{node}=inv(p(i).C(ix,ix));
    end

    if i < levels
        p(i).omega=1.0;    % smoother should not overrelax
        p(i).M = interp_m(p(i).m,p(i).n,ratio,ones(2)); % mask for coarsening the obstacles
        p(i+1).m=(p(i).m-1)/ratio+1; % prepare next level
        p(i+1).n=(p(i).n-1)/ratio+1;
        p(i+1).dx=p(i).dx*ratio;
        p(i+1).dy=p(i).dy*ratio;
        p(i+1).bs=bs;
        p(i+1).nnods=p(i+1).m*p(i+1).n;
        p(i+1).ndofs=p(i+1).bs*p(i+1).nnods;
        p(i+1).zlim=[-(p(i).zlim(2)-p(i).zlim(1))/2,(p(i).zlim(2)-p(i).zlim(1))/2]; % display limits
        p(i+1).penalty=p(i).penalty*4;
    else
        p(i).omega=1.0;  % with constraints, omega>1 may cause divergence
    end
end


for mg_it=1:20,
    for i=1:levels-1
        % pre-smoothing
        disp_solution(i,p(i),sprintf('level %g it %g start',i,mg_it))
        for sm_it=1:5    % need at least one iteration to hopefully satisfy obstacles
            [p(i).a,p(i).err]=relax_constr(p(i).a,p(i));
            disp_solution(i,p(i),sprintf('level %g mg it %g pre-smoothing it %g',i,mg_it,sm_it))
        end
        % coarsening
        r = p(i).b - p(i).C * p(i).a;  % residual
        p(i+1).b = p(i+1).Ir.' * r;  % coarsen residual
        p(i+1).l= spmat_trans_maxvec(p(i).M,p(i).l - p(i).a(1:bs:end));  % coarsen lower bound, <=0
        p(i+1).u=-spmat_trans_maxvec(p(i).M,p(i).a(1:bs:end) - p(i).u);  % coarsen upper bound, >=0
        p(i+1).a = zeros(p(i+1).ndofs,1);             % initial approx
    end
    % coarsest solve
    i=levels;
    
    p(i).err=[];
    disp_solution(i,p(i),sprintf('level %g it %g start',i,mg_it))
    
    for sm_it=1:100
        [p(i).a,p(i).err]=relax_constr(p(i).a,p(i));
        disp_solution(i,p(i),sprintf('level %g mg it %g coarsest solve it %g',i,mg_it,sm_it))
        if p(i).err.diff <= norm(p(i).a,inf)*1e-12,
            break,
        end
    end
    for i=levels-1:-1:1
        % interpolate the correction
        p(i).a = p(i).a + p(i+1).I*p(i+1).a;       
        p(i).err.err=[];
        disp_solution(i,p(i),sprintf('level %g mg it %g after coarse correction',i,mg_it))
 
        % post-smoothing
        for sm_it=1:10
            [p(i).a,p(i).err]=relax_constr(p(i).a,p(i));
            disp_solution(i,p(i),sprintf('level %g mg it %g post-smoothing in %g',i,mg_it,sm_it))
        end
    end
end
p(1).a_disp=reshape(p(1).a(1:bs:end),[m,n]);
