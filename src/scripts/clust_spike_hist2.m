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

% subsets = {'subset1' 'subset2'} ; 
subsets = {'subset1'} ; 


for sdx = subsets

    spk_avghist.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
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

indat1 = spk_avghist.mean.subset1(1:20,cortmask)' ;
indat2 = spk_avghist.mean.subset2(1:20,cortmask)';

nodesz = sqrt(2* size(indat1,1)+0.25)-0.5+1 ; 

%%

[uu,vv,~,~,v] = pca(diff(indat1')','NumComponents',4) ; 

%%


cc = corr(indat1) ; 

%imagesc(indat1(sortedInd(vv(:,4)),:))

N = 20 ; 
qq = nan(N,1) ; 
parti = ones(N,1) ;
parti(1) = 2 ; 
for idx = 1:N-1
    [~,qq(idx)] = fcn_get_q(cc,parti,'potts') ; 
    parti(idx+1) = 2 ; 
end
