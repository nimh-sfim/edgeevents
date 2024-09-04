function [ edgenodeinds ] = get_ets_node_inds(nodevec) 

nodevec = nodevec(:) ; 

% if sum(nodevec)~=1
%     error('not good input')
% end

nnodes = length(nodevec) ;

% get the upper triangle mask
trium = logical(triu(ones(nnodes),1)) ; 
% and the inds of the maskffff
triui = find(trium) ; 

inds = find((nodevec*~nodevec')+(nodevec*~nodevec')') ; 
edgenodeinds = ismember(triui,inds) ; 

end