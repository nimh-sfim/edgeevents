%% start

clearvars
clc
close all

%% preamble load data

%run('./config/config_hcp_sch200_1.m') 

scannames = {'REST1_LR' 'REST1_RL'} ; 

fmridat = struct() ; 
for idx = scannames
    fmridat.(idx{1}) = load_hcp_alldata(...
        [ DD.INTERM '/hcp352_nusregts_FIX2phys_schaefer200/' ],...
        'schaefer200-yeo17',...
        sublist.all, ...
        ['*' idx{1} '*']) ; 
end

%% make variablity maps 

% subsets = {'subset1' 'subset2'} ; 

dvarsdat  = struct() ; 
kurtdat = struct() ;
stddat  = struct() ; 
meandat = struct() ; 
tsnrdat = struct() ; 

for sdx = scannames

    dvarsdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ; 
    kurtdat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ; 
    stddat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ; 
    meandat.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ;  

    tsnrmap.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ; 
    acmap.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.all)) ; 

    for idx = 1:length(sublist.all)
    
        disp(idx)
    
        %sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        tmp = fmridat.(sdx{1})(idx).ts(:,1:finfo.nnodes) ;

        dvarsdat.(sdx{1}).map(:,idx) = rms(diff(tmp)) ;
        kurtdat.(sdx{1}).map(:,idx) = kurtosis(tmp) ; 
        stddat.(sdx{1}).map(:,idx) = std(tmp) ; 
        meandat.(sdx{1}).map(:,idx) = mean(tmp) ; 
        tsnrdat.(sdx{1}).map(:,idx) = mean(tmp)./std(tmp) ; 
        
        acmap.(sdx{1}).map(:,idx) = get_acf_hwhm(zscore(tmp),0.72,15) ;  
    end

end

%% simply save this data for future use... 

save('./data/interim/ts_variability.mat',"kurtdat","dvarsdat","stddat","meandat","tsnrdat","acmap")

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

% %% compare to spk conn
% 
% filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
% load(filename)
% 
% %%
% sss = mean(stdmats.subset1.mats,3) ; 
% imsc_grid_comm(sss,parc.ca(1:200))
% axis square
% colorbar
% colormap(inferno(100))
% 
% %%
% 
% clf
% tiledlayout(1,4)
% 
% names = {'short' 'inter' 'long'}
% 
% for idx = 1:3
%     nexttile()
%     scatter_w_rho(tv(spike_conn.subset1.(names{idx})),tv(sss))
%     title([ names{idx} ' counts vs edgestd'])
%     ylabel('std x std')
%     xlabel([ names{idx} ' counts'])
% end
% 
% nexttile()
% scatter_w_rho(tv(meanfc),tv(sss))
% title([ 'meanfc vs edgestd'])


