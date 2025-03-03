function [blockmat,phigh,plow,permblockdat] = run_blocky_spintest(indat,ca,perms,nperm) 

blockmat = get_blocky(indat,ca,@(x_) mean(x_,'all','omitmissing')) ; 
permsz = size(perms,2) ; 

nblocks = size(blockmat,1) ; 
permblockdat = nan(nblocks,nblocks,nperm) ; 
for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perms(:,randi(permsz)) ; 
    permblockdat(:,:,idx) = get_blocky(indat(ii,ii),ca,@(x_) mean(x_,'all','omitmissing')) ; 
end

permdathigh = sum(blockmat<=permblockdat,3) ; 
permdatlow = sum(blockmat>=permblockdat,3) ; 

phigh = (permdathigh+1)./(nperm+1) ; 
plow = (permdatlow+1)./(nperm+1) ;
