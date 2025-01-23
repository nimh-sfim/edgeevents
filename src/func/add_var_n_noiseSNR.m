function outts = add_var_n_noiseSNR(ints,varfrac,SNR)

if nargin < 3
    error('need three args')
end

maxamp = max(abs(ints)) ;

% add variability 
tt1 = generate_phase_surrogates(ints); 
tt1 = tt1./(max(abs(tt1))).*varfrac.*maxamp ; 

% add gaus noise
tt2 = (mean(ints + tt1)/SNR).*randn(size(ints)) ; 

outts = ints + tt1 + tt2 ; 