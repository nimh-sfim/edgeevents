%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

% %% hubs from full dataset
% 
% [meanfc_mask,ii] = my_omst(1-meanfc) ; 
% cortexplot(sum(meanfc_mask))

% %% get "functional hub" from the meanfc matrix
% 
% ts = datStr(15).ts(:,1:200) ; 
% ets = get_ets(ts) ; 
% fc = corr(ts) ; 
% [~,sortfc] = sort(tv(fc),'descend') ; 

%%

% prct = 95 ; 
prcts = [80 95 99] ;
funcconsit = zeros(length(sublist.subset1),finfo.nnodes,length(prcts)) ; 
funcdegree = zeros(length(sublist.subset1),finfo.nnodes) ; 

for idx = 1:length(sublist.subset1)

    disp(idx)

    ts = datStr(idx).ts(:,1:finfo.nnodes) ; 
    ets = get_ets(ts) ; 
    fc = corr(ts) ; 

    for pdx = 1:length(prcts)
        prct = prcts(pdx) ; 
    
        tsabvthr = ts>prctile(ts,prct) ; 
        etsabvthr = ets_2_ts(ets>prctile(ets,prct)) ; 
    
        tt = tsabvthr .* etsabvthr ; 

        %funcleaderscr = sum(etsabvthr.*tsabvthr) ; 
        funcconsit(idx,:,pdx) = arrayfun(@(i_) mean(nonzeros((etsabvthr(:,i_).*(tsabvthr(:,i_))))),1:size(ts,2)) ; 
    end 

    % tmp = fc>prctile(fc,80) ; 
    % tmp = tmp | tmp' ; 

    cc = corr(ts) ; 
    funcdegree(idx,:) = mean(cc.*(cc>0)) ; 

end

%%

% group gradient
% addpath('~/joshstuff/git_pull/BrainSpace/matlab/')
% gg = GradientMaps('kernel','na','approach', 'Laplacian Eigenmap') ; 
% ggg = gg.fit(meanfc) ; 

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

%%

parc_plot_wcolorbar(mean(funcconsit(:,:,1)),surfss,annotm,...
    [min(mean(funcconsit(:,:,1))) max(mean(funcconsit(:,:,1)))],parula(100),[100 100 600 1000])

parc_plot_wcolorbar(mean(funcdegree),surfss,annotm,...
    [min(mean(funcdegree)) max(mean(funcdegree))],parula(100),[100 100 600 1000])


%% 

scatter_w_rho(mean(funcdegree),mean(funcconsit(:,:,1)))
xlabel('functional degree')
ylabel('functional spk consistency')
axis square

%%

scatter_w_rho(grad_cifti(1,:),mean(funcconsit(:,:,1)))
xlabel('functional grad 1')
ylabel('functional spk consistency')
axis square

%% functional consistency given canonical systems???

