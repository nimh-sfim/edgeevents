%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

high_bin = 4 ; 
maxspk = 1200  ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

spklen_names = {'short' 'inter' 'long'} ; 

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

% fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

refparc = parc.ca(1:200) ; 

nedges = ( finfo.nnodes * (finfo.nnodes-1) ) / 2 ; 

%% now read in spike lengths to make a histogram of lengths

for sdx = subsets

    spike_rss.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    overall_rss.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        ets = get_ets(datStr(sind).ts) ; 
        ets = ets(:,cortmask) ; 
            
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        % origfc = corr(datStr(sind).ts(:,1:200)) ; 

        %% lets get the 3 rss 
        tmprss = zeros(finfo.ntp,3) ; 
        for hdx = 1:3
            tmprss(:,hdx) = get_rss(dd==hdx) ;
            

        end

        spike_rss.(sdx{1}){idx} = tmprss ; 
        overall_rss.(sdx{1}){idx} = get_rss(readdat.spike_len_mat(:,cortmask)) ; 

    end
end

%% 

rrr = load_hcp_regressors('/Users/faskowitzji/joshstuff/data/hcp352_regressors/hcp352_regressors',...
      'rfMRI_REST1_RL','*rfMRI_REST1_RL_confounds.tsv') ; 

%% 

% res1 = zeros(101,length(sublist.(sdx{1}))) ; 
% res2 = zeros(101,length(sublist.(sdx{1}))) ; 
% res3 = zeros(101,length(sublist.(sdx{1}))) ; 

cres1 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
cres2 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
cres3 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
cc = zeros(3,1) ; 


ooo1 = zeros(1196,19900) ; 
ooo2 = zeros(1196,19900) ; 
ooo3 = zeros(1196,19900) ; 

for idx = 1:length(sublist.(sdx{1}))

    disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 

    sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 

    ts = datStr(sind).ts ; 
    ets = get_ets(ts) ; 
    ets = ets(:,cortmask) ; 
    
    %sig = rrr(datStr(sind).sub).values{1}.framewise_displacement ; 
    % sig = mean(ts,2) ; 

    filename = [DD.PROC '/' 'REST1_RL' '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
    readdat = load(filename,'spike_len_mat') ; 
    
    d = readdat.spike_len_mat(:,cortmask) ; 

    % simvec = 1:10 ; 
    % 
    % for sdx = 1:simvec
    % 
    % end

    dd = discretize( d, lowmedhigh_edges) ; 
    dd(isnan(dd)) = 0 ;
    dd = dd(3:end-2,:) ; 

    [~,o1] = count_spks(dd==1) ; 
    [~,o2] = count_spks(dd==2) ; 
    [~,o3] = count_spks(dd==3) ; 

    ooo1 = ooo1 + o1 ; 
    ooo2 = ooo2 + o2 ; 
    ooo3 = ooo3 + o3 ; 

    % [~,o1spk] = findpeaks(rms(o1,2),...
    %     "MinPeakProminence",std(rms(o1,2))*2,...
    %     "MinPeakHeight",prctile(rms(o1,2),80)) ; 
    % o1spk = ismember(1:finfo.ntp,o1spk) ; 
    % 
    % [~,o2spk] = findpeaks(rms(o2,2),...
    %     "MinPeakProminence",std(rms(o2,2))*2,...
    %     "MinPeakHeight",prctile(rms(o2,2),80)) ; 
    % o2spk = ismember(1:finfo.ntp,o2spk) ; 
    % 
    % [~,o3spk] = findpeaks(rms(o3,2),...
    %     "MinPeakProminence",std(rms(o3,2))*2,...
    %     "MinPeakHeight",prctile(rms(o3,2),80)) ; 
    % o3spk = ismember(1:finfo.ntp,o3spk) ; 
    % 
    % allspk = [ o1spk; o2spk; o3spk ]' ; 
    % 
    % ee = get_entdrops([get_rss((o1)) ...
    %     get_rss((o2)) get_rss((o3))]) ; 
    % 
    % % ee = get_entdrogps(...
    % %     normalize([get_rss((o1)) get_rss((o2)) get_rss((o3))],'range',[0 1])...
    % %     ,20) ; 
    % 
    % % rsscols = [rms(o1,2) rms(o2,2) rms(o3,2)] ; 
    % % 
    % % ee = get_entdrops(...
    % %     rsscols ...
    % %     ,20) ; 
    % 
    % % [ee,~] = get_entdrops(...
    % %     normalize([get_rss((dd==1)) get_rss((dd==2)) get_rss((dd==3))],'norm',2),...
    % %     20) ; 
    % 
    % % cc = cc + countcount(ee,1:3) ; 
    % % 
    % % cres1(:,:) = sum(cat(3,cres1(:,:),mksq(mean(ets(ee==1,:),1))),3,'omitnan') ; 
    % % cres2(:,:) = sum(cat(3,cres2(:,:),mksq(mean(ets(ee==2,:),1))),3,'omitnan') ; 
    % % cres3(:,:) = sum(cat(3,cres3(:,:),mksq(mean(ets(ee==3,:),1))),3,'omitnan') ; 
    % 
    % 
    % cres1(:,:,idx) = cov(ets_2_ts(o1)) ;  
    % cres2(:,:,idx) = cov(ets_2_ts(o2)) ;  
    % cres3(:,:,idx) = cov(ets_2_ts(o3)) ;  

end

%%

node2sysload = @(ts_,ca_) cell2mat(arrayfun(@(i_) mean(ts_(:,ca_==i_),2), 1:max(ca_),'UniformOutput',false)) ; 

nn = node2sysload(ets_2_ts(o3),parc.ca(1:200))' ; 
make_trans_dens(nn,10)
[~,mi] = max(nn) ; 
mi(sum(nn)==0) = 0 ; 



%% TODOOOOOO
% look at the cov(ets_2_ts(o2) of the spike onset times for the different
% spike lengths... maybe there's different connectivity at the different
% spike onsets...
% and maybe run an identifiability analysis or something



%%

r1 = overall_rss.subset1{1} ; 
r2 = spike_rss.subset1{1} ; 

diffpeaks = ...
    findpeaks(zscore(r2(:,3))-zscore(r2(:,1)), ...
    'MinPeakHeight',prctile(zscore(r2(:,3))-zscore(r2(:,1)),80), ...
    'MinPeakProminence',std(zscore(r2(:,3))-zscore(r2(:,1)))*2.5) ;


s1 = get_ets_states(readdat.spike_len_mat(:,cortmask)>0,parc.ca(1:200),'assortweak') ; 
tmpets = get_ets(datStr(sind).ts) ; 
s2 = get_ets_states(tmpets(:,cortmask),parc.ca(1:200),'assortweak') ; 

%%  

c1 = corr(ts,get_rss(dd==1)) ; 
c2 = corr(ts,get_rss(dd==2)) ; 
c3 = corr(ts,get_rss(dd==3)) ; 
aa = get_acf_hwhm(zscore(ts),0.72) ; 

e1 = arrayfun(@(i_)myent_simp(d(i_,:),1:20),1:finfo.ntp) ; 
e2 = arrayfun(@(i_)myent_simp(dd(i_,:),1:20),1:finfo.ntp) ; 

st1 = arrayfun(@(i_)std(d(i_,:)),1:finfo.ntp) ; 
st2 = arrayfun(@(i_)std(dd(i_,:)),1:finfo.ntp) ; 


etsrss = get_rss(ets) ; 
spkrss = get_rss()

%%

ts = datStr(sind).ts ; 
ets = get_ets(ts) ; 
ets = ets(:,cortmask) ; 

%sig = rrr(datStr(sind).sub).values{1}.framewise_displacement ; 
% sig = mean(ts,2) ; 

filename = [DD.PROC '/' 'REST1_LR' '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
readdat = load(filename,'spike_len_mat') ; 

d = readdat.spike_len_mat(:,cortmask) ; 
dd = discretize( d, lowmedhigh_edges) ; 
dd(isnan(dd)) = 0 ;
% 
[~,o1] = count_spks(dd==1) ; 
[~,o2] = count_spks(dd==2) ; 
[~,o3] = count_spks(dd==3) ; 

nn = normalize([ get_rss(ets_2_ts(o1)) get_rss(ets_2_ts(o2)) get_rss(ets_2_ts(o3)) ],2,'norm',1) ;
entts = (-sum(nn.*log2(nn),2,'omitnan')) ; 
entdrops = entts<prctile(entts,20) ; 
[ ~, entst ] = max(nn(entdrops,:),[],2) ;
entdrops_annot = zeros(length(entdrops),1) ; 
entdrops_annot(entdrops) = entst ; 