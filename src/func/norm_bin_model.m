function [zout,xx,mm] = norm_bin_model(refdata,refbins,indat)
% zout is the zscore values within the bins
% the xx and mm are to be used for plotting the point in the middle of bin

nbins = length(refbins)-1 ; 

refdata_disc = discretize(refdata,refbins) ; 

mm = zeros(nbins,1)  ; 

% xx = movmean(refbins(2:end),[1 0]) ; 
% xx = [ xx(2)-(xx(3)-xx(2)) refbins(2:(end-1)) ] ; % xx(end-1)-xx(end-2)+xx(end-1) ] ; 

refbins2 = refbins.*1 ; 
refbins2(1) = min(refdata(:)) ; 
refbins2(end) = max(refdata(:)) ;

xx = interp1(1:(nbins+1),refbins2,1.5:(nbins+0.5)) ;

% xx = refbins(2:(end-1)) + diff(refbins(2:(end)))
% xx(end) = max(refdata(:)) - xx(end-1) ; 
% xx = []

zout = zeros(size(refdata)) ; 
% for each of the bins, zscore
for idx = 1:nbins
    
    ee = refdata_disc==idx ; 
    [zz,mm(idx)] = zscore(indat(ee)) ; 
    zout(ee) = zz ; 

end