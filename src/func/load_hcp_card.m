function [outmap,outmap2] = load_hcp_card(dataDir,sname) 

if nargin < 2
    error('need 2 args')
end

%%


dd = dir([dataDir '/*' sname '*tsv' ]) ; 

outmap = dictionary() ; 
outmap2 = dictionary() ; 
tmpscan = sname ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop it
for idx = 1:length(dd) 
    
    tmpdat = table2array(readtable([ dd(idx).folder '/' dd(idx).name ],...
        'FileType','text','Delimiter','\t')) ; 
    tmpsub = regexprep(dd(idx).name,'_card.*','') ; 

    tmpsig = resampsig1d(tmpdat(:,1),25,1/0.72) ; 

    if (idx == 1) || ~isKey(outmap,tmpsub) 
        outmap(tmpsub) = dictionary(tmpscan,{tmpsig}) ; 
        outmap2(tmpsub) = dictionary(tmpscan,{tmpdat(:,1)}) ; 
    else
        tt = outmap(tmpsub) ; 
        tt(tmpscan) = {tmpsig}  ;
        outmap(tmpsub) = tt ; 

        tt = outmap2(tmpsub) ; 
        tt(tmpscan) = {tmpdat(:,1)}  ;
        outmap2(tmpsub) = tt ; 
    end

end