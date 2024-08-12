function scount = count_spks(indat,thr)

if nargin < 2
    spks = indat ; 
else
    spks = indat>thr ;
end

scount = sum(...
    zerocrossrate(...
    double(spks),... % needs to be double here!!! 
    'Level',1,'WindowLength',1,'TransitionEdge','rising')) ; 

