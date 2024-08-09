function [lenmat,spklens,spkinds] = spk_lenmat(spk)
% replace spikes with their length values 

[spklens,lencell,spkinds] = arrayfun(@(i_) get_contact_times(spk(:,i_)) , 1:size(spk,2) , ...
    'UniformOutput',false) ; 
lenmat = cell2mat(lencell) ; 
