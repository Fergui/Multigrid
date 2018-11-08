function mg=multigrid_method(u,p)
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
levels=levels-1
[m,n]=size(u);
% initialize level 1
mg(1).sm=sm;
mg(1).sn=sn;
[mg(1).Sm,mg(1).Sn]=meshgrid(sm,sn);
mg(1).m=m;
mg(1).n=n;
mg(1).x=1:mg(1).m;
mg(1).y=1:mg(1).n;
mg(1).xc=1:ratio:mg(1).m;
mg(1).yc=1:ratio:mg(1).n;
mg(1).dx=p.dx;
mg(1).dy=p.dy;
mg(1).nnods=mg(1).m*mg(1).n;
[Hs,gs]=interpolation(mg(1).Sm,mg(1).Sn,p.shapes,times);
Hc=cat(1,Hs{:});
[H,rows]=condence(Hc);
gc=cat(1,gs{:});
g=gc(rows);
mg(1).H=H;
mg(1).g=g;
R=p.R(sm,sn);
mg(1).R=R;
mg(1).q=p.q;
mg(1).nmesh=p.nmesh;
mg(1).max_step=p.max_step;
mg(1).max_depth=p.max_depth;
mg(1).min_depth=p.min_depth;
mg(1).select=p.select;
mg(1).ofunc=p.ofunc;
mg(1).vmask=p.vmask(sm,sn);
mg(1).penalty=p.penalty;
% initial fire arrival time
mg(1).ui=u;
% initial contribution matrix
[mg(1).Ci,mg(1).Ji]=vJ(mg(1).ui,mg(1).R,mg(1));

% prepare all the static things in each level
for i=2:levels
    mg(i).sm=mg(i-1).sm(1:ratio:mg(i-1).m); 
    mg(i).sn=mg(i-1).sn(1:ratio:mg(i-1).n);
    [mg(i).Sm,mg(i).Sn]=meshgrid(mg(i).sm,mg(i).sn);
    mg(i).m=length(mg(i).sm);
    mg(i).n=length(mg(i).sn);
    mg(i).x=mg(i-1).xc;
    mg(i).y=mg(i-1).yc;
    mg(i).xc=1:ratio:mg(i).m;
    mg(i).yc=1:ratio:mg(i).n;
    mg(i).dx=mg(i-1).dx*ratio;
    mg(i).dy=mg(i-1).dy*ratio;
    mg(i).nnods=mg(i).m*mg(i).n;
    [Hs,gs]=interpolation(mg(i).Sm,mg(i).Sn,p.shapes,times);
    Hc=cat(1,Hs{:});
    [mg(i).H,rows]=condence(Hc);
    gc=cat(1,gs{:});
    mg(i).g=gc(rows);
    mg(i).R=mg(i-1).R(mg(i-1).xc,mg(i-1).yc);
    mg(i).q=p.q;
    mg(i).nmesh=p.nmesh;
    mg(i).max_step=p.max_step;
    mg(i).max_depth=p.max_depth;
    mg(i).min_depth=p.min_depth;
    mg(i).select=p.select;
    mg(i).ofunc=p.ofunc;
    mg(i).vmask=mg(i-1).vmask(mg(i-1).xc,mg(i-1).yc);
    mg(i).penalty=p.penalty;
    mg(i).ui=mg(i-1).ui(mg(i-1).xc,mg(i-1).yc);
    [mg(i).Ci,mg(i).Ji]=vJ(mg(i).ui,mg(i).R,mg(i));
end

exp=0:levels-1;
strategy=2.^exp;

for i=levels:-1:1
    figure(i),
    if i<levels,
        mg(i).u=interp2(mg(i+1).Sm,mg(i+1).Sn,mg(i+1).u,mg(i).Sm,mg(i).Sn);
        [mg(i).C,mg(i).Jop]=vJ(mg(i).u,mg(i).R,mg(i));
        Js=mg(i).Jop;
    else
        mg(i).u=mg(i).ui;
        mg(i).Jop=mg(i).Ji;
        Js=mg(i).Jop;
    end
    subplot(2,2,1), mesh(mg(i).u), view([0,1]), title(['Level ',num2str(i),': First approximation']), drawnow;
    subplot(2,2,2), pcolor(mg(i).R), title(['Level ',num2str(i),': Rate of Spread']), drawnow;
    for sm_it=1:strategy(i)
        fprintf('Level %d, iteration %d.\n',i,sm_it);
        tic
        dir=zeros(mg(i).m,mg(i).n);
        for ii=1:mg(i).m
            for jj=1:mg(i).n
                if mg(i).vmask(ii,jj)==1
                    dir(ii,jj)=1;
                    for sig=-1:2:1
                        d=sig*dir;
                        [Jmin,~,um,~]=linesearch(mg(i).u,@gJ,mg(i).R,d,mg(i));
                        if Jmin<mg(i).Jop(end)
                            mg(i).Jop=[mg(i).Jop;Jmin];
                            mg(i).u=um;
                        end
                    end
                    dir(ii,jj)=0;
                end
            end
        end
        Js=[Js,mg(i).Jop(end)];
        subplot(2,2,3), plot(Js), title(['Level ',num2str(i),': Objective function value']), drawnow;
        subplot(2,2,4), mesh(mg(i).u), view([0,1]), title(['Level ',num2str(i),': Solution at iteration ',num2str(sm_it)]), drawnow;
        toc
    end
end
end
