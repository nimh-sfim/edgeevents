function [gmean,phigh,plow,permmeandat] = run_groupmean_spintest(indat,ca,perms,nperm) 

permsz = size(perms,2) ; 
ngroup = length(unique(ca)) ; 
%n = size(indat,1) ; 
uca = unique(ca(:)) ; 

gmeanfunc = @(dat_,ca_) arrayfun(@(i_) mean(dat_(ca_==i_),'omitnan'),uca) ; 
% can get away with uca here cause ca dones't change

gmean = gmeanfunc(mean(indat,2,'omitnan'),ca) ; 

permmeandat = nan(ngroup,nperm) ; 
for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perms(:,randi(permsz)) ; 

    tmp = indat(ii,ii) ; 
    permmeandat(:,idx) = gmeanfunc(mean(tmp,2,'omitnan'),ca) ; 
end

permdathigh = sum(gmean<=permmeandat,2) ; 
permdatlow = sum(gmean>=permmeandat,2) ; 

phigh = (permdathigh+1)./(nperm+1) ; 
plow = (permdatlow+1)./(nperm+1) ;
