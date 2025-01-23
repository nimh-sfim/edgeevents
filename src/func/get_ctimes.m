function [meanIET,maxIET,stdIET] = get_ctimes(inSpk) 
% assume in the inSpk is ntpXchannels, and that spikes are 1's 

% get cell array of contact times
ctimes = arrayfun(@(i_) get_contact_times(~inSpk(:,i_)), 1:size(inSpk,2), ...
    'UniformOutput',false) ; 
meanIET = cellfun(@(i_) mean(i_) , ctimes ) ; 
maxIET = cellfun(@(i_) max(i_) , ctimes ) ; 
stdIET = cellfun(@(i_) std(i_) , ctimes ) ; 
