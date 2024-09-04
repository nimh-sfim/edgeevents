function [Pmat,Pnorm,state_vals] = make_trans_prob(symbol_seq,state_sz)

if nargin < 2
    state_sz = 'uniq' ;
end

ntp = length(symbol_seq) ;

switch state_sz
    case 'max'
        state_sz = max(symbol_seq) ;
        state_vals = 1:state_sz ; 
        rectsymbseq = symbol_seq ; 
    case 'uniq'
        state_vals = unique(symbol_seq) ; 
        state_sz = length(state_vals) ;

        rectsymbseq = nan(ntp,1) ;
        for idx = 1:state_sz
            tmp = state_vals(idx) ;
            rectsymbseq(symbol_seq==tmp) = idx ;
        end
    otherwise % the user has specified a number
        state_vals = 1:state_sz ;
        rectsymbseq = symbol_seq ; 
end

Pmat = zeros(state_sz) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% populate matrix
s1 = rectsymbseq(1) ;
for idx = 1:(ntp-1)
    s2 = rectsymbseq(idx+1) ;
    Pmat(s1,s2) = Pmat(s1,s2) + 1;
    s1 = s2 ;
end

% normalize rows to add to one!
Pnorm = normalize(Pmat,2,"norm",1) ;
