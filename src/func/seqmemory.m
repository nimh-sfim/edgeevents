function [ mem ] = seqmemory(sequence,corrtype) 
% one lag correlation
% n.b. biased estimator

if nargin < 2
    corrtype = 'p' ;
end

if ~isvector(sequence) 
    error('input needs to be vector')
end

if isscalar(sequence)
    mem = nan ;
    return 
end

n = length(sequence) ;
mem = corr(sequence(1:(n-1)),sequence(2:n),'type',corrtype) ; 
