function [ out ] = ets_2_ts(ets,func)

if nargin < 2
    func = 'mean' ; 
end

gathermethod = 1 ; 
switch func
    case 'mean'
        funchandle = @(x_) mean(x_,2) ; 
    case 'sum'
        funchandle = @(x_) sum(x_,2) ; 
    case 'ec'
        funchandle = @(x_) eigenvector_centrality_und(mks(x_)) ;
        gathermethod = 2 ;
    case 'cc'
        funchandle = @(x_) clustering_coef_wu(mks(x_)) ; 
        gathermethod = 2 ;
    otherwise
        error('not an option for func')
end


if gathermethod == 1 
    ne = size(ets,2) ; 
    nn = sqrt(2 * ne + 0.25) - 0.5 + 1 ; 
    
    out = cell2mat(...
        arrayfun(@(i_) ...
        funchandle(ets(:,get_ets_node_inds((1:nn)==i_))), ...
        1:nn, 'UniformOutput',false ) ...
        ) ; 

elseif gathermethod == 2
    out = cell2mat(arrayfun(@(i_) ... 
        funchandle(ets(i_,:)), 1:size(ets,1), ...
        'UniformOutput',false))' ; 

end

end

function [ a ] = mks(M)
% make triu + diag outta vect

K = 1 ;

b = length(M) ;
s = sqrt(2 * b + 0.25) - 0.5 + K;
a = triu(ones(s),K) ;
a(logical(a)) = M ;
a = a + triu(a,K)' ; 

end 
