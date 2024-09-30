function [outmap] = load_hcp_regressors(dataDir,sname,fname) 

if nargin < 3
    error('need 3 args')
end

%%


dd = dir([dataDir '/*/*' sname '*/' fname ]) ; 

outmap = dictionary() ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop it
for idx = 1:length(dd) 

    [a,tmpscan] = fileparts(dd(idx).folder) ; 
    [~,tmpsub] = fileparts(a) ; 
    tmpdat = readtable([ dd(idx).folder '/' dd(idx).name ],...
        'FileType','text','Delimiter','\t') ; 
    
    if (idx == 1) || ~isKey(outmap,tmpsub) 
        outmap(tmpsub) = dictionary(tmpscan,{tmpdat}) ; 
    else
        tt = outmap(tmpsub) ; 
        tt(tmpscan) = {tmpdat}  ;
        outmap(tmpsub) = tt ; 
    end

end