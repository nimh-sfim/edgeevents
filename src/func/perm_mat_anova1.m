function [anovastuff,tbl,stats] = perm_mat_anova1(mat,labs,nperms)

if nargin < 3
    nperms = 5000 ;
end

N = size(mat,1) ; 
trium = logical(triu(ones(N),1)) ; 
[u,v] = find(trium) ; 
ca_u = labs(u) ; 
ca_v = labs(v) ; 

ulabs = unique(labs) ; 
% nuniq = length(ulabs) ; 

% inline
getmatinds = @(i_,i_ca_u,i_ca_v) logical(triu(mksq(logical((i_ca_u == i_) & (i_ca_v == i_))),1)) ; 

% run real anova
yy = cell2mat(arrayfun(@(i_) mat(getmatinds(i_,ca_u,ca_v)),ulabs,'UniformOutput',false)) ; 
gg = cell2mat(arrayfun(@(i_) ones(sum(getmatinds(i_,ca_u,ca_v),'all'),1).*i_,ulabs,'UniformOutput',false)) ; 
[~,tbl,stats] = anova1(yy,gg,'off') ; 
anovastuff.f = tbl{2,5} ; 
anovastuff.df_report = [tbl{2,3} tbl{3,3}] ; 

% now the permuation to get p on F val
pfvals = nan(nperms,1); 
for idx = 1:nperms
    
    plabs = labs(randperm(length(labs))) ; 
    p_ca_u = plabs(u) ; 
    p_ca_v = plabs(v) ; 

    pyy = cell2mat(arrayfun(@(i_) mat(getmatinds(i_,p_ca_u,p_ca_v)),ulabs,'UniformOutput',false)) ; 
    pgg = cell2mat(arrayfun(@(i_) ones(sum(getmatinds(i_,p_ca_u,p_ca_v),'all'),1).*i_,ulabs,'UniformOutput',false)) ;
    [~,ptbl,~] = anova1(pyy,pgg,'off') ; 
    pfvals(idx) = ptbl{2,5} ; 

end % perm loop

anovastuff.p = (sum(pfvals>=anovastuff.f)+1)/(nperms+1) ; 

end

function [ a , asym ] = mksq(M,K)
% make triu + diag outta vect

if nargin < 2 
   K = 1 ; 
end

b = length(M) ;
s = sqrt(2 * b + 0.25) - 0.5 + K;
a = triu(ones(s),K) ;
a(logical(a)) = M ;
a = a + triu(a,K)' ; 

% make the asymmetric too
if nargout > 1
    asym = triu(a,K) ; 
end

end