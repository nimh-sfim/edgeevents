function [ent] = myent_simp(x_,c_)

if nargin < 2
    c_ = unique(x_) ; 
end

co = countcount(x_(:),c_(:)) ;
p = co ./ sum(co) ;
lp = log2(p) ;
lp(isinf(lp)) = 0 ; 

ent = -1 * sum(p.*lp) ; 
if isnan(ent)
    ent = 0 ; 
end