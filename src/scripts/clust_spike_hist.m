%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

spk_bins = [1:20 100] ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 


for sdx = subsets

    spk_avghist.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize(readdat.spike_len_mat,spk_bins)  ;
        dd(isnan(dd)) = 0 ; 
    
        ss = cell2mat(arrayfun(@(i_) count_spks(dd==i_),1:length(spk_bins),UniformOutput=false)') ;
        spk_avghist.(sdx{1}){idx} = ss ;  

        % get a spk dist for each edge 
        % ss = cell2mat(arrayfun(@(i_) countcount(dd(:,i_),1:21),1:size(dd,2),'UniformOutput',false)) ; 
      
    end

end

%% lets make average

for sdx = subsets

    spk_avghist.mean.(sdx{1}) = zeros(length(spk_avghist),255*254/2) ; 
    spk_avghist.nmean1.(sdx{1}) = zeros(length(spk_avghist),255*254/2) ; 
    spk_avghist.nmean2.(sdx{1}) = zeros(length(spk_avghist),255*254/2) ; 

    for idx = 1:length(sublist.(sdx{1}))

        disp(idx)

        tmp = spk_avghist.mean.(sdx{1}) ; 
        spk_avghist.mean.(sdx{1}) = tmp + spk_avghist.(sdx{1}){idx} ; 

        tmp = spk_avghist.nmean1.(sdx{1}) ; 
        tmp2 = normalize(spk_avghist.(sdx{1}){idx},'norm',1) ; 
        tmp2(isnan(tmp2)) = 0 ; 
        spk_avghist.nmean1.(sdx{1}) = tmp + tmp2 ;  


        tmp = spk_avghist.nmean2.(sdx{1}) ; 
        tmp2 = normalize(spk_avghist.(sdx{1}){idx},'norm',2) ; 
        tmp2(isnan(tmp2)) = 0 ; 
        spk_avghist.nmean2.(sdx{1}) = tmp + tmp2 ;  

        tmp = spk_avghist.nmean2.(sdx{1}) ; 
        tmp2 = normalize(spk_avghist.(sdx{1}){idx},'norm',2) ; 
        tmp2(isnan(tmp2)) = 0 ; 
        spk_avghist.nmean2.(sdx{1}) = tmp + tmp2 ;  


    end

    spk_avghist.mean.(sdx{1}) = spk_avghist.mean.(sdx{1}) ./ length(sublist.(sdx{1})) ; 
    spk_avghist.nmean1.(sdx{1}) = spk_avghist.nmean1.(sdx{1}) ./ length(sublist.(sdx{1})) ;
    spk_avghist.nmean2.(sdx{1}) = spk_avghist.nmean2.(sdx{1}) ./ length(sublist.(sdx{1})) ;

end

%% save it

filename = [ DD.PROC '/spk_hist_' OUTSTR '.mat' ] ; 
save(filename,'spk_avghist','spk_bins','-v7.3')

%% run the kmeans in a consensus manner

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;
% mask again to exlcude limbic 
nolimb = find(ismember(parc.ca(1:finfo.nnodes),find(cellfun(@(x_)~isempty(x_),regexp(parc.names,'Limb'))))) ; 
nolimbmask = ~ismember(u,nolimb) & ~ismember(v,nolimb) ; 

nodesnolimb = ~ismember(1:200,nolimb) ; 
namesnolimb = parc.ca(1:200) ; 
namesnolimb = namesnolimb(~ismember(namesnolimb,find(cellfun(@(x_)~isempty(x_),regexp(parc.names,'Limb'))))) ; 

%% initial sweep

sweep_iters = 100 ; 
sweep_k = 2:12 ; 

indat1 = spk_avghist.nmean1.subset1(1:20,cortmask)' ;
indat2 = spk_avghist.nmean1.subset2(1:20,cortmask)';

nodesz = sqrt(2* size(indat1,1)+0.25)-0.5+1 ; 

% mykm = @(k_)(kmeans(indat1, k_,...
%     'Replicates',5,'Distance','sqeuclidean'));
mykm = @(k_)(kmedoids(indat1, k_,'Algorithm','clara',...
    'Replicates',25,'Distance','cityblock'));

sweep_iter_ca = zeros(nodesz*(nodesz-1)/2,length(sweep_k),sweep_iters) ; 

eval_k = struct() ; 
eval_names = {'CalinskiHarabasz' 'DaviesBouldin' 'silhouette' } ; 
for jdx = eval_names 
    eval_k.subset1.(jdx{1}) = cell(sweep_iters,1) ; 
    eval_k.subset2.(jdx{1}) = cell(sweep_iters,1) ; 
end

for idx = 1:sweep_iters

    tt = tic() ; 

    % make clusterings over sweep_k 
    sweep_iter_ca(:,:,idx) = cell2mat(arrayfun(@(i_) mykm(i_) , sweep_k  , 'UniformOutput' , false)) ;  

    for jdx = eval_names 
        eval_k.subset1.(jdx{1}){idx} = evalclusters(indat1,sweep_iter_ca(:,:,idx),jdx{1}) ; 
        eval_k.subset2.(jdx{1}){idx} = evalclusters(indat2,sweep_iter_ca(:,:,idx),jdx{1}) ; 
    end

    elapsed_time = toc(tt) ; 
    disp([ 'iter ', num2str(idx) ,' : ' num2str(round(elapsed_time,2)) ' sec' ])

end

% aggregrate via versatility

sweep_vers = zeros(nodesz*(nodesz-1)/2,length(sweep_k)) ; 

for idx = 1:length(sweep_k)
    disp(idx)

    sweep_vers(:,idx) = node_versatility(squeeze(sweep_iter_ca(:,idx,:))) ; 
end

% find centroid
sweep_k_vi = zeros(sweep_iters,sweep_iters,length(sweep_k)) ; 
sweep_k_visum = zeros(sweep_iters,length(sweep_k),1) ; 
% sweep_k_sil = zeros(sweep_iters,length(sweep_k),1) ; 

for idx = 1:length(sweep_k)
    disp(idx)

    sweep_k_vi(:,:,idx) = partition_distance(squeeze(sweep_iter_ca(:,idx,:))) ; 
    sweep_k_visum(:,idx) = sum(sweep_k_vi(:,:,idx)) ;
    % sweep_k_sil(:,idx) = arrayfun(@(i_) mean(silhouette(indat1,sweep_iter_ca(:,idx,i_))),1:sweep_iters) ; 
end

% find the centroid at each resolution 

sweep_maxsil_ind = zeros(length(sweep_k),1) ; 
sweep_minvi_ind = zeros(length(sweep_k),1) ;

for idx = 1:length(sweep_k)

    [~,sweep_maxsil_ind(idx)] = min(sweep_k_sil(:,idx)) ;
    [~,sweep_minvi_ind(idx)] = min(sum(sweep_k_vi(:,:,idx) )) ;

end

% save the data we need for further analysis 
filename = [ DD.PROC '/spike_hist_clust.mat' ] ;
save(filename,'sweep_*', '-v7.3')


% %% lets try running pca??
% 
% [uu,vv,~,~,ee] = pca(normalize(spk_avghist.mean.subset1(1:20,cortmask&nodesnolimb)',2,'norm',2),...
%     'NumComponents',5) ; 
% 
% [uu,vv,~,~,ee] = pca(spk_avghist.mean.subset1(1:20,cortmask)',...
%     'NumComponents',5) ; 
% 
% 
% %%
% 
% ii = spk_avghist.mean.subset1(1:20,cortmask)' ; 
% ii2 = normalize(normalize(spk_avghist.mean.subset1(1:20,cortmask)','norm',2),2,'norm',1) ; 
% 
% for idx = repmat(1:20,1,100) 
%     tiledlayout(1,2)
%     nexttile()
%     imsc_grid_comm(mksq(ii(:,idx)),parc.ca(1:200)) ; title(idx) ; 
%     nexttile()
%     imsc_grid_comm(mksq(ii2(:,idx)),parc.ca(1:200)) ; title(idx) ; 
%     waitforbuttonpress
% end
% 
% %%  
% 
% for idx = repmat(1:20,1,100) 
% 
%     imsc_grid_comm(mksq(spk_avghist.mean.subset1(idx,cortmask)'),parc.ca(1:200)) ; title(idx) ; colorbar
%     clim([0 3])
%     waitforbuttonpress
% end
% 
% %%
% 
% ecs = zeros(200,20) ; 
% for idx = 1:20
%     ecs(:,idx) = eigenvector_centrality_und(mksq(spk_avghist.mean.subset1(idx,cortmask)')) ; 
% end
% 
% 
% cortexplot(mean(ecs(:,4:20),2))
% 
% cortexplot(mean(ecs(:,1:2),2))

