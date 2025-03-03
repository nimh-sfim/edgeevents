function [phigh,plow] = run_mat_spintest(indat,perms,nperm) 

permsz = size(perms,2) ; 
n = size(indat,1) ; 

permdathigh = zeros(n) ; 
permdatlow = zeros(n) ; 

for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perms(:,randi(permsz)) ; 

    tmp = indat(ii,ii) ; 
    tmp(1:n+1:end) = nan ; 

    permdathigh = permdathigh + (indat<=tmp) ; 
    permdatlow = permdatlow + (indat>=tmp) ; 

end

phigh = (permdathigh+1)./(nperm+1) ; 
plow = (permdatlow+1)./(nperm+1) ;
