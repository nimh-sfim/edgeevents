function [out] = get_rss(in,meanit)
if nargin < 2
    meanit = false ;
end

if ~meanit
    out = sqrt(sum(in.^2,2)) ; 
else
    out = sqrt(mean(in.^2,2)) ; 
end
