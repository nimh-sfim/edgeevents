function [burstVec,memVec,meanIET,maxIET] = spk_burst_mem_2(inSpk) 
% assume in the inSpk is ntpXchannels, and that spikes are 1's 

% get cell array of contact times
ctimes = arrayfun(@(i_) get_contact_times(~inSpk(:,i_)), 1:size(inSpk,2), ...
    'UniformOutput',false) ; 

% burstiness corrected for finite sequences
bform2 = @(y) (sqrt(length(y)+1)*(std(y)/mean(y)) - sqrt(length(y)-1))/((sqrt(length(y)+1)-2)*(std(y)/mean(y)) + sqrt(length(y)-1));

burstVec = cellfun(@(i_) bform2(i_) , ctimes ) ; 
burstVec(burstVec==-1) = NaN ; % address times where there's only one spike

memVec = cellfun(@(i_) seqmemory(i_) , ctimes ) ; 

if nargout > 2
    meanIET = cellfun(@(i_) mean(i_) , ctimes ) ; 
end

if nargout > 3
    maxIET = cellfun(@(i_) max(i_) , ctimes ) ; 
end
