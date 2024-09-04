%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%% make a tsnr map

subsets = {'subset1' 'subset2'} ; 

dvarsdat  = struct() ; 
kurtdat = struct() ;
vardat  = struct() ; 
tsnrdat = struct() ; 

for sdx = subsets

    dvarsdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    kurtdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    vardat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    tsnrdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        tmp = datStr(sind).ts(:,1:finfo.nnodes) ;

        dvarsdat.(sdx{1}).map(:,idx) = rms(diff(tmp)) ;
        kurtdat.(sdx{1}).map(:,idx) = kurtosis(tmp) ; 
        vardat.(sdx{1}).map(:,idx) = var(tmp) ; 
        tsnrdat.(sdx{1}).map(:,idx) = mean(tmp) ; 

    end

end

%% simply save this data for future use... 

save('./data/interim/ts_variability.mat',"kurtdat","dvarsdat","vardat")
