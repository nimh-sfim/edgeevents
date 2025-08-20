function [out] = get_simple_hrv_rmssd(cardsig,hz_signal,wsz)

ntp = length(cardsig) ; 
t_vec = [0:1/hz_signal:(1/hz_signal*(ntp-1))]; %current time signal


[~,pksind] = findpeaks(cardsig,'MinPeakProminence',std(cardsig),...
    'MinPeakDistance',hz_signal/1.8) ; % max 108 bpm, hack

tmp = zeros(length(cardsig),1) ; 
tmp(pksind(1:end-1)) = diff(t_vec(pksind)) ; 

winInds = (0:1:ntp-1)' + (1:wsz) ;
% adjust to center the window, lazy way?
winInds = winInds - floor(wsz/2) ; 

% calculate rows that exceed the ntp, or are less than 0
winInds(winInds>ntp) = 0 ;
winInds(winInds<1) = 0 ;

out = nan(ntp,1) ; 
for idx = 1:size(winInds,1)
    iii = winInds(idx,:) ; 
    out(idx) = rmssd(tmp(nonzeros(iii))) ; 
end