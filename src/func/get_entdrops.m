function [entdrops,entts] = get_entdrops(indat,droppct)

if nargin < 2
    droppct = 20 ; 
end

indatnorm = normalize(indat,2,'norm',1) ; % normalize along 1st dim, sum of 1
entts = (-sum(indatnorm.*log2(indatnorm),2,'omitnan')) ; % entropy calc
% entts(isnan(indatnorm)) = nan ; 
entts(isnan(sum(indatnorm,2))) = nan ; 
entdrops_vec = entts<=prctile(entts,droppct) ; % get the drops, will exclude nans
[ ~, entst ] = max(indatnorm(entdrops_vec,:),[],2) ; % get colum with max val
entdrops = zeros(length(entdrops_vec),1) ; 
entdrops(entdrops_vec) = entst ; % get the columns at each entdrop moment 
