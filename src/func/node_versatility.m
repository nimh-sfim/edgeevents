function [V] = node_versatility(ca)
% ca should be nodes x paritions

CM = agreement(ca) ./ size(ca,2);

Cs = sin(pi*CM);
V = sum(Cs, 1)./sum(CM, 1);
V(V<1e-10) = 0;
