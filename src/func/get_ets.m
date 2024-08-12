function [ets,etsten] = get_ets(ts) 

ets = fcn_ets(ts) ;
if nargout > 1
    etsten = get_etsten(ts,0) ;
end

end

function [a,aa] = fcn_ets(ts)

[~,n] = size(ts);               % number samples/nodes
z = zscore(ts);                 % z-score
[u,v] = find(triu(ones(n),1));  % get edges
a = z(:,u).*z(:,v);             % edge ts products
if nargout>1
    aa = sqrt(z(:,u).^2 + z(:,v).^2); 
end
end

function etsten = get_etsten(ts,diagval)

if nargin < 2
    diagval = 1 ;
end

n = size(ts,2) ;
z = zscore(ts);                 % z-score
etsten = permute(...
    bsxfun(@times,permute(z,[1,2,3]),permute(z,[1,3,2])),...
    [ 2 3 1 ]) ;

if isnan(diagval)
    return
else % set the diagonal values to the diagval
    etsten((repmat(eye(n),[ 1 1 size(etsten,3)]))==1) = diagval ; 
end


end % end main
