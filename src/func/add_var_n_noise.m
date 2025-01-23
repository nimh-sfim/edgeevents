function outts = add_var_n_noise(ints,varfrac,noisefrac)

if nargin < 3
    error('need three args')
end

maxamp = max(abs(ints)) ;

% add variability 
tt1 = generate_phase_surrogates(ints); 
tt1 = tt1./(max(abs(tt1))).*varfrac.*maxamp ; 

% add uniform noise
tt2 = (rand(size(ints))-0.5).*2.*noisefrac.*maxamp ;

outts = ints + tt1 + tt2 ; 