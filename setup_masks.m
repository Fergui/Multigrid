function p=setup_masks(p)
%p=setup_masks(p)
% Set up the masks mask and vmask in the structure p
%input:
%   p    structure with:
%           per1_mask   matrix, mask of the first perimeter
%           per2_mask   matrix, mask of the second perimeter
%output:
%   p    structure with:
%           mask        matrix, true where the values of level set function have to change
%           vmask       matrix, true where the values of level set function can change
[m,n]=size(p.per2_mask);
% mask = 1 where the level set function is computed, 0 where given
p.mask=true(m,n);
p.mask(p.per1_mask)=0;
p.mask(p.per2_mask)=0;
% vmask = 1 where the equation is evaluated: add per2 boundary
p.vmask=p.mask;
neigh=[0 0; 0 1; 1 0; 0 -1; -1 0];
nneigh=size(neigh,1);
[m,n]=size(p.mask);
for i=1:m
    for j=1:n
        if p.per2_mask(i,j),
            for k=1:nneigh
                ii=i+neigh(k,1);
                jj=j+neigh(k,2);
                if ii>=1 || i<=m || jj>=1 || jj<=n || ~p.mask(ii,jj),
                    p.vmask(i,j)=0;
                    break
                end
            end
        end
    end
end
end