function [scount,spkedge] = count_spks(indat,thr)

if nargin < 2
    spks = indat ; 
else
    spks = indat>thr ;
end

spkedge = zerocrossrate(...
    double(spks),... % needs to be double here!!! 
    'Level',1,'WindowLength',1,'TransitionEdge','rising')>0 ;

scount = sum(spkedge) ; 
