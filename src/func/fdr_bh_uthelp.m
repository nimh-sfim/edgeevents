function outpass_mat = fdr_bh_uthelp(pmat,varargin) 

n = size(pmat,1) ; 
utmask = logical(triu(ones(n))) ; 
outpass = fdr_bh(pmat(utmask),varargin{:}) ; 
outpass_mat = zeros(n) ; 
outpass_mat(utmask) = outpass ; 
