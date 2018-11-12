function mg=multigrid_method_opt(u,p)
if ~isfield(p,'tol')
    p.tol=1e-6;
end
if ~isfield(p,'coarse')
    p.coarse=0;
end
if ~isfield(p,'dynR')
    p.dynR=0;
end
times=unique(p.g);
[mp,np]=size(u);
sm=findmgfit(mp);
sn=findmgfit(np);
u=u(sm,sn);
% setup
ratio=2;  
[mm,nn]=size(u);
% determine number of levels
for levels=1:100,
    mm=(mm-1)/ratio+1; 
    nn=(nn-1)/ratio+1;
    if mm~=round(mm) || nn~=round(nn) || mm<4 || nn<4
        break
    end
end
levels=3
[m,n]=size(u);

%% Initialize level 1: finest level
% Dimensions and resolutions
mg(1).m=m;
mg(1).n=n;
mg(1).nnods=mg(1).m*mg(1).n;
mg(1).dx=p.dx;
mg(1).dy=p.dy;
% Grid and indexes
mg(1).ii=1:mg(1).m;
mg(1).jj=1:mg(1).n;
mg(1).iic=1:ratio:mg(1).m;
mg(1).jjc=1:ratio:mg(1).n;
mg(1).X=p.X(sm,sn);
mg(1).Y=p.Y(sm,sn);
mg(1).Xc=mg(1).X(mg(1).iic,mg(1).jjc);
mg(1).Yc=mg(1).Y(mg(1).iic,mg(1).jjc);
% H and g for the linear system of constraints
[Hs,gs]=interpolation(mg(1).X,mg(1).Y,p.shapes,times);
Hc=cat(1,Hs{:});
[H,rows]=condence(Hc);
gc=cat(1,gs{:});
g=gc(rows);
mg(1).H=H;
mg(1).g=g;
% Rate of spread
R=p.R(sm,sn);
mg(1).R=R;
% Other configuration parameters
mg(1).q=p.q;
mg(1).nmesh=p.nmesh;
mg(1).max_step=p.max_step;
mg(1).max_depth=p.max_depth;
mg(1).min_depth=p.min_depth;
mg(1).select=p.select;
mg(1).ofunc=p.ofunc;
mg(1).vmask=p.vmask(sm,sn);
mg(1).penalty=p.penalty;
% Initial fire arrival time
mg(1).ui=u;
% Initial contribution matrix
[mg(1).Ci,mg(1).Ji]=vJ(mg(1).ui,mg(1).R,mg(1));
% ROS dynamic important variables
if p.dynR
    ffuel=fields(p.nfuelcat);
    for k=1:length(ffuel)
        mg(1).nfuelcat.(ffuel{k})=p.nfuelcat.(ffuel{k})(sm,sn);
    end
    mg(1).ignS.dzdxf=p.ignS.dzdxf(sm,sn);
    mg(1).ignS.dzdyf=p.ignS.dzdyf(sm,sn);
    mg(1).dynS.simt=p.dynS.simt;
    mg(1).dynS.uf=p.dynS.uf(sm,sn,:);
    mg(1).dynS.vf=p.dynS.vf(sm,sn,:);
    mg(1).dynS.fmc_g=p.dynS.fmc_g(sm,sn,:);
end

%% Prepare all the static things in each level
for i=2:levels
    % Dimensions and resolutions
    mg(i).m=length(mg(i-1).iic);
    mg(i).n=length(mg(i-1).jjc);
    mg(i).dx=mg(i-1).dx*ratio;
    mg(i).dy=mg(i-1).dy*ratio;
    mg(i).nnods=mg(i).m*mg(i).n;
    % Grid and indexes
    mg(i).ii=mg(i-1).iic; 
    mg(i).jj=mg(i-1).jjc;
    mg(i).iic=1:ratio:mg(i).m;
    mg(i).jjc=1:ratio:mg(i).n;
    mg(i).X=mg(i-1).Xc;
    mg(i).Y=mg(i-1).Yc;
    mg(i).Xc=mg(i).X(mg(i).iic,mg(i).jjc);
    mg(i).Yc=mg(i).Y(mg(i).iic,mg(i).jjc);
    % H and g for the linear system of constraints
    [Hs,gs]=interpolation(mg(i).X,mg(i).Y,p.shapes,times);
    Hc=cat(1,Hs{:});
    [mg(i).H,rows]=condence(Hc);
    gc=cat(1,gs{:});
    mg(i).g=gc(rows);
    % Rate of spread
    if p.coarse
        mg(i).R=coarsening(mg(i-1).R,2);
        mg(i).R(~mg(i).vmask)=0;
    else
        mg(i).R=mg(i-1).R(mg(i).ii,mg(i).jj);
    end
    % Other configuration parameters
    mg(i).q=p.q;
    mg(i).nmesh=p.nmesh;
    mg(i).max_step=p.max_step;
    mg(i).max_depth=p.max_depth;
    mg(i).min_depth=p.min_depth;
    mg(i).select=p.select;
    mg(i).ofunc=p.ofunc;
    mg(i).vmask=mg(i-1).vmask(mg(i).ii,mg(i).jj);
    %mg(i).vmask=(coarsening(mg(i-1).vmask,2)>0);
    mg(i).penalty=p.penalty;
    % Initial fire arrival time
    if p.coarse
        mg(i).ui=coarsening(mg(i-1).ui,2);
    else
        mg(i).ui=mg(i-1).ui(mg(i).ii,mg(i).jj);
    end
    % Initial contribution matrix
    [mg(i).Ci,mg(i).Ji]=vJ(mg(i).ui,mg(i).R,mg(i));
    % ROS dynamic important variables
    if p.dynR
        ffuel=fields(p.nfuelcat);
        for k=1:length(ffuel)
            mg(i).nfuelcat.(ffuel{k})=mg(i-1).nfuelcat.(ffuel{k})(mg(i).ii,mg(i).jj);
        end
        mg(i).ignS.dzdxf=mg(i-1).ignS.dzdxf(mg(i).ii,mg(i).jj);
        mg(i).ignS.dzdyf=mg(i-1).ignS.dzdyf(mg(i).ii,mg(i).jj);
        mg(i).dynS.simt=mg(i-1).dynS.simt;
        mg(i).dynS.uf=mg(i-1).dynS.uf(mg(i).ii,mg(i).jj,:);
        mg(i).dynS.vf=mg(i-1).dynS.vf(mg(i).ii,mg(i).jj,:);
        mg(i).dynS.fmc_g=mg(i-1).dynS.fmc_g(mg(i).ii,mg(i).jj,:);
    end
end

%% Strategy vector
exp=0:levels-1;
strategy=2*2.^exp;
strategy(end)=500;

%% Cascadic Multigrid
for i=levels:-1:1
    if i<levels,
        if p.dynR
            interpolantR=scatteredInterpolant(mg(i+1).X(:),mg(i+1).Y(:),mg(i+1).R(:));
            %interpolantR.Method='natural';
            Rn=interpolantR(mg(i).X(:),mg(i).Y(:));
            mg(i).R=reshape(Rn,mg(i).m,mg(i).n);
            mg(i).R(~mg(i).vmask)=0;
        end
        interpolantu=scatteredInterpolant(mg(i+1).X(:),mg(i+1).Y(:),mg(i+1).u(:));
        %interpolantu.Method='natural';
        un=interpolantu(mg(i).X(:),mg(i).Y(:));
        mg(i).u=reshape(un,mg(i).m,mg(i).n);
        [mg(i).C,mg(i).Jop]=vJ(mg(i).u,mg(i).R,mg(i));
        Js=mg(i).Jop;
    else
        mg(i).u=mg(i).ui;
        mg(i).C=mg(i).Ci;
        mg(i).Jop=mg(i).Ji;
        Js=mg(i).Jop;
        if p.dynR
            mg(i).ds=dyninterp(mg(i).u,mg(i));
            mg(i).R=ros_file(mg(i).u,mg(i).ignS,mg(i).ds,mg(i));
            mg(i).R(~mg(i).vmask)=0;
        end
    end
    U=mg(i).u; U(~mg(i).vmask)=nan;
    figure(i), subplot(2,3,1), mesh(U), view([0,1]), title(['Level ',num2str(i),': First approximation']), drawnow;
    figure(i), subplot(2,3,2), mesh(mg(i).R), view([0,1]), title(['Level ',num2str(i),': Rate of Spread']), colorbar, drawnow;
    figure(i), subplot(2,3,3), h=pcolor(mg(i).vmask); title(['Level ',num2str(i),': Mask']), set(h,'EdgeColor','None'), drawnow;
    ri=1:mg(i).m;
    rj=1:mg(i).n;
    dir=zeros(mg(i).m,mg(i).n);
    for sm_it=1:strategy(i)
        figure(i), subplot(2,3,4), h=pcolor(mg(i).C); title(['Level ',num2str(i),': Contribution matrix']), set(h,'EdgeColor','None'), colorbar, drawnow;
        fprintf('Level %d, iteration %d.\n',i,sm_it);
        tic
        for ii=ri
            for jj=rj
                if mg(i).vmask(ii,jj)==1
                    li=max(1,ii-2); ui=min(ii+2,mg(i).m);
                    lj=max(1,jj-2); uj=min(jj+2,mg(i).n);
                    il=li:ui;
                    jl=lj:uj;
                    dir(ii,jj)=1;
                    for sig=-1:2:1
                        d=projectdir(sig*dir,mg(i).H);
                        [Jmin,~,um,Jlow,Cm,pm]=linesearch_multigrid(mg(i).u(il,jl),@vJ,mg(i).R(il,jl),d(il,jl),mg(i).vmask(il,jl),mg(i));
                        if Jmin<Jlow
                            C=mg(i).C; 
                            C(li:ui-2,lj:uj-2)=Cm;
                            Jn=(sum(C(:).^mg(i).q)*mg(i).dx*mg(i).dy)^(1/mg(i).q)+pm;
                            if Jn<mg(i).Jop(end)
                                mg(i).C=C;
                                mg(i).Jop=[mg(i).Jop;Jn];
                                mg(i).u(il,jl)=um;
                            end
                        end
                     end
                    dir(ii,jj)=0;
                end
            end
        end
        Js=[Js,mg(i).Jop(end)];
        figure(i), subplot(2,3,5), plot(Js,'-*'), title(['Level ',num2str(i),': Objective function value']), drawnow;
        U=mg(i).u; U(~mg(i).vmask)=nan;
        figure(i), subplot(2,3,6), mesh(U), view([0,1]), title(['Level ',num2str(i),': Solution at iteration ',num2str(sm_it)]), drawnow;
        toc
        if Js(end-1)-Js(end)<eps
            break
        end
        ri=flip(ri);
        rj=flip(rj);
        if p.dynR
            mg(i).ds=dyninterp(mg(i).u,mg(i));
            mg(i).R=ros_file(mg(i).u,mg(i).ignS,mg(i).ds,mg(i));
            mg(i).R(~mg(i).vmask)=0;
            figure(i), subplot(2,3,2), mesh(mg(i).R), view([0,1]), title(['Level ',num2str(i),': Rate of Spread']), colorbar, drawnow;
        end
    end
end
UI=mg(1).ui; UI(~mg(1).vmask)=nan;
U=mg(1).u; U(~mg(1).vmask)=nan;
figure, 
subplot(1,2,1), mesh(UI), view([0,1]), title('Initial approximation');
subplot(1,2,2), mesh(U), view([0,1]), title('Final result after the cascadic multigrid');
end