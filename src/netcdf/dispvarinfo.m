function dispvarinfo(p)
% dispvarinfo(p)
% display info on one variable
% 
% Jan Mandel
%-------------------------------------------------------------------------

    if isfield(p,'var_value'),
        disp([p.varname,' size ',num2str(size(p.var_value)),' netcdf stored as ',...
              p.vartype_m,' (',num2str(p.dimlength),')'])
    end
end
