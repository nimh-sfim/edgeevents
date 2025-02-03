function [burst1,burst2,memVec,meanIET,maxIET] = spk_burst_mem(inSpk) 
% assume in the inSpk is ntpXchannels, and that spikes are 1's 
%
% burstiness: -1=regular, 0=poisson, >0=bursty
%   definition corrected for finite sequences is unbounded positive

% get cell array of contact times
ctimes = arrayfun(@(i_) get_contact_times(~inSpk(:,i_)), 1:size(inSpk,2), ...
    'UniformOutput',false) ; 

% burstiness corrected for finite sequences
bform2 = @(y) (sqrt(length(y)+1)*(std(y)/mean(y)) - sqrt(length(y)-1))/((sqrt(length(y)+1)-2)*(std(y)/mean(y)) + sqrt(length(y)-1));
% orig burstiness
bform1 = @(y) (std(y)-mean(y)) / (std(y)+mean(y));

burst2 = cellfun(@(i_) bform2(i_) , ctimes ) ; 
burst2(burst2==-1) = NaN ; % address times where there's only one spike

% burstiness corrected for finite sequences

burst1 = cellfun(@(i_) bform1(i_) , ctimes ) ; 
burst1(burst1==-1) = NaN ; % address times where there's only one spike

memVec = cellfun(@(i_) seqmemory(i_) , ctimes ) ; 

if nargout > 3
    meanIET = cellfun(@(i_) mean(i_) , ctimes ) ; 
end

if nargout > 4
    maxIET = cellfun(@(i_) max(i_) , ctimes ) ; 
end
