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
stddat = struct() ; 
meandat = struct() ; 

for sdx = subsets

    dvarsdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    kurtdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    vardat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    stddat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    meandat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        tmp = datStr(sind).ts(:,1:finfo.nnodes) ;

        dvarsdat.(sdx{1}).map(:,idx) = rms(diff(tmp)) ;
        kurtdat.(sdx{1}).map(:,idx) = kurtosis(tmp) ; 
        vardat.(sdx{1}).map(:,idx) = var(tmp) ; 
        stddat.(sdx{1}).map(:,idx) = std(tmp) ; 
        meandat.(sdx{1}).map(:,idx) = mean(tmp) ; 

    end

end

%% simply save this data for future use... 

save('./data/interim/ts_variability.mat',"kurtdat","dvarsdat","vardat")

%%

parc_plot_wcolorbar(mean(stddat.(sdx{1}).map,2),surfss,annotm,...
    [0 max(mean(stddat.(sdx{1}).map,2))],inferno(100),[100 100 600 1000])


%% make std mult maps

stdmats = struct() ; 

for sdx = subsets

    stdmats.(sdx{1}).mats = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))

        ss = stddat.(sdx{1}).map(:,idx) ; 
        stdmats.(sdx{1}).mats(:,:,idx) = ss*ss' ; 

    end
end

%% compare to spk conn
 
filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
load(filename)

%%
sss = mean(stdmats.subset1.mats,3) ; 
imsc_grid_comm(sss,parc.ca(1:200))
axis square
colorbar
colormap(inferno(100))

%%

clf
tiledlayout(1,4)

names = {'short' 'inter' 'long'}

for idx = 1:3
    nexttile()
    scatter_w_rho(tv(spike_conn.subset1.(names{idx})),tv(sss))
    title([ names{idx} ' counts vs edgestd'])
    ylabel('std x std')
    xlabel([ names{idx} ' counts'])
end

nexttile()
scatter_w_rho(tv(meanfc),tv(sss))
title([ 'meanfc vs edgestd'])



