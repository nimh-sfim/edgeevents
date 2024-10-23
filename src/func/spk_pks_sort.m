function [pks,locs] = spk_pks_sort(dat,varargin)

if ~isvector(dat)
    error('needs vector')
end


[pks,locs] = findpeaks(dat,varargin{:}) ;
% sort the peaks by mag
[~,ss] = sort(pks,'descend') ; 
pks = pks(ss) ; 
locs = locs(ss) ; 

