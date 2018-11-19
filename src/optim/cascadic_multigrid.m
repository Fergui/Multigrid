function mg=cascadic_multigrid(mg,strategy)
levels=length(mg);
%% Cascadic Multigrid
for i=levels:-1:1
    if i<levels,
        if mg(i).dynR
            if mg(i).interp=='cubic'
                mg(i).R=interp2(mg(i+1).II,mg(i+1).JJ,mg(i+1).R,mg(i).II,mg(i).JJ,'cubic');
            elseif mg(i).interp=='akima'
                mg(i).R=interp2(mg(i+1).II,mg(i+1).JJ,mg(i+1).R,mg(i).II,mg(i).JJ,'akima');
            elseif mg(i).interp=='linear'
                interpolantR=scatteredInterpolant(mg(i+1).X(:),mg(i+1).Y(:),mg(i+1).R(:));
                Rn=interpolantR(mg(i).X(:),mg(i).Y(:));
                mg(i).R=reshape(Rn,mg(i).m,mg(i).n);
            else
                error('Error: Any valid interpolation method assigned!');
            end
            mg(i).R(~mg(i).vmask)=0;
        end
        if mg(i).interp=='cubic'
            mg(i).u=interp2(mg(i+1).II,mg(i+1).JJ,mg(i+1).u,mg(i).II,mg(i).JJ,'cubic');
        elseif mg(i).interp=='akima'
            mg(i).u=interp2(mg(i+1).II,mg(i+1).JJ,mg(i+1).u,mg(i).II,mg(i).JJ,'akima');
        elseif mg(i).interp=='linear'
            interpolantu=scatteredInterpolant(mg(i+1).X(:),mg(i+1).Y(:),mg(i+1).u(:));
            un=interpolantu(mg(i).X(:),mg(i).Y(:));
            mg(i).u=reshape(un,mg(i).m,mg(i).n);
        else
            error('Error: Any valid interpolation method assigned!');
        end
        [mg(i).C,mg(i).Jop]=vJ(mg(i).u,mg(i).R,mg(i));
        Js=mg(i).Jop;
    else
        if mg(i).dynR
            mg(i).ds=dyninterp(mg(i).u,mg(i));
            mg(i).R=ros_file(mg(i).u,mg(i).ignS,mg(i).ds,mg(i));
            mg(i).R(~mg(i).vmask)=0;
        end
        [mg(i).C,mg(i).Jop]=vJ(mg(i).u,mg(i).R,mg(i));
        Js=mg(i).Jop;
    end
    U=mg(i).u; U(~mg(i).vmask)=nan;
    figure(i), subplot(2,3,1), mesh(U), view([0,1]), title(['Level ',num2str(i),': First approximation']), drawnow;
    figure(i), subplot(2,3,2), mesh(mg(i).R), view([0,1]), title(['Level ',num2str(i),': Rate of Spread']), colorbar, drawnow;
    figure(i), subplot(2,3,3), h=pcolor(mg(i).vmask); title(['Level ',num2str(i),': Mask']), set(h,'EdgeColor','None'), drawnow;
    ri=1:mg(i).m;
    rj=1:mg(i).n;
    for sm_it=1:strategy(i)
        fprintf('Level %d, iteration %d.\n',i,sm_it);
        figure(i), subplot(2,3,4), h=pcolor(mg(i).C); caxis([-3,3]); title(['Level ',num2str(i),': Contribution matrix']), set(h,'EdgeColor','None'), colorbar, drawnow;
        tic
        dir=zeros(mg(i).m,mg(i).n);
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
        if mg(i).dynR
            mg(i).ds=dyninterp(mg(i).u,mg(i));
            mg(i).R=ros_file(mg(i).u,mg(i).ignS,mg(i).ds,mg(i));
            mg(i).R(~mg(i).vmask)=0;
            figure(i), subplot(2,3,2), mesh(mg(i).R), view([0,1]), title(['Level ',num2str(i),': Rate of Spread']), colorbar, drawnow;
        end
    end
end
end