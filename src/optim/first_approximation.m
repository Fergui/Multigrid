function u=first_approximation(p,objf)
[p.C,p.Jop]=objf(p.u,p.R,p);
Js=p.Jop;
U=p.u; U(~p.vmask)=nan;
figure(100), subplot(2,3,1), mesh(U), view([0 1]),  title(['First approximation']), drawnow;
for k=1:p.max_iter
    if k==1
        M=zeros(p.m,p.n);
    end
    figure(100), subplot(2,3,2), h=pcolor(p.C); caxis([0,1]), title(['Contribution matrix']), set(h,'EdgeColor','None'), colorbar, drawnow;
    dir=zeros(p.m,p.n);
    for ii=1:p.m
        for jj=1:p.n
            if p.vmask(ii,jj)==1
                li=max(1,ii-2); ui=min(ii+2,p.m);
                lj=max(1,jj-2); uj=min(jj+2,p.n);
                il=li:ui;
                jl=lj:uj;
                dir(ii,jj)=1;
                dd=max(projectdir(dir,p.H),0);
                if k==1
                    M=M+dd;
                end
                for sig=-1:2:1
                    d=sig*dd;
                    [Jmin,~,um,Jlow,Cm,pm]=linesearch_multigrid(p.u(il,jl),objf,p.R(il,jl),d(il,jl),p.vmask(il,jl),p);
                    if Jmin<Jlow
                        C=p.C; 
                        C(li:ui-2,lj:uj-2)=Cm;
                        Jn=sum(C(:))+pm;
                        if Jn<p.Jop(end)
                            p.C=C;
                            p.Jop=[p.Jop;Jn];
                            p.u(il,jl)=um;
                        end
                    end
                end
                dir(ii,jj)=0;
            end
        end
    end
    figure(100), subplot(2,3,3), plot(p.Jop,'-*'),  title(['Objective function after each linesearch']), drawnow;
    if k==1
        figure(100), subplot(2,3,4), h=pcolor(M); title(['Points modified']), set(h,'EdgeColor','None'), colorbar, drawnow;
    end
    Js=[Js,p.Jop(end)];
    figure(100), subplot(2,3,5), plot(Js,'-*'),  title(['Objective function after each iteration']), drawnow;
    U=p.u; U(~p.vmask)=nan;
    figure(100), subplot(2,3,6), mesh(U), view([0 1]),  title(['New solution']), drawnow;
end
u=p.u;

end