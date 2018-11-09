function [H,g] = interpolation(X,Y,shapes,times)
    ns=length(shapes);
    H=cell(ns,1);
    g=cell(ns,1);
    for s=1:length(shapes)
        Hs=interop_bary(X,Y,shapes(s).x,shapes(s).y);
        [H{s},~]=condence(Hs);
        g{s}=times(s)*ones(size(H{s},1),1);
    end
end