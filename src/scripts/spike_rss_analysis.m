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

    spike_rsspeaks.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 


    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        % origfc = corr(datStr(sind).ts(:,1:200)) ; 

        %% lets get the 3 rss 
        tmprss = zeros(finfo.ntp,3) ; 
        tmppeaks = zeros(finfo.ntp,3) ; 
        for hdx = 1:3
            tmprss(:,hdx) = get_rss(dd==hdx) ;
            [~,b] = findpeaks(tmprss(:,hdx),'MinPeakHeight',prctile(tmprss(:,hdx),80),...
                'MinPeakProminence',std(tmprss(:,hdx))*2) ; 
            tmppeaks(:,hdx) = ismember(1:finfo.ntp,b) ; 
        end

        spike_rsspeaks.(sdx{1}){idx} = tmppeaks ;  

    end
end

%% now go thru spikes and get pattern corresponding to spike frames

for sdx = subsets

    spkpttrn_sizes.(sdx{1}) = sum(cell2mat(spike_rsspeaks.(sdx{1}))) ; 
    
    spkpttrns.(sdx{1}) = cell(3,1) ; 
    
    % initialize 
    for pdx = 1:3
        spkpttrns.(sdx{1}){pdx} = zeros(spkpttrn_sizes.(sdx{1})(pdx),nedges) ; 
    end
    
    % populate the big mats
    counter = ones(3,1) ; 

    for idx = 1:length(sublist.(sdx{1}))
   
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
   
        tmpets = get_ets(datStr(sind).ts) ; 

        for pdx = 1:3

            % get the spikes for this subject, at this len
            spkvec = spike_rsspeaks.(sdx{1}){idx}(:,pdx) ; 
            tmpsz = sum(spkvec) ; 
            
            % the patterns
            spkpttrns.(sdx{1}){pdx}(counter(pdx):(counter(pdx)+tmpsz-1),:) = ...
                single(tmpets(logical(spkvec),cortmask)) ; 

            counter(pdx) = counter(pdx) + tmpsz ; 
        end

    end

end

%% save the data

% filename = [ DD.PROC '/spk_rsspeaks_' OUTSTR '.mat' ] ; 
% save(filename,'spkpttrns','spike_rsspeaks')

%% now make the big comparisons

for sdx = subsets
    
    spkpttrncmp.(sdx{1}) = cell(3,1) ; 
    
    for pdx = 1:3
  
        tt = tic ; 
        disp(['making comparison ' num2str(pdx) ' ' sdx{1}])  ; 

        linmat = single(lins_ccc(spkpttrns.(sdx{1}){pdx}')) ; 

        % take the condorance between all pairs of patterns
        spkpttrncmp.(sdx{1}){pdx} = linmat ; 
        toc(tt)

    end
end

%% save the lins mats

filename = [ DD.PROC '/spk_rsspeakssim_' OUTSTR '.mat' ] ; 
save(filename,'spkpttrncmp')

%% do clustering of these 

mod_reps = 500 ; 

for sdx = subsets
    
    spkpttrnclust.(sdx{1}).mods = cell(3,1) ; 
    
    for pdx = 1:3
  
        clusterdat = spkpttrncmp.(sdx{1}){pdx} ; 

        mods = nan(size(clusterdat,1),mod_reps) ; 
        ss = sdx{1} ; 

        for idx = 1:mod_reps
            disp(['mod rep' num2str(idx) ' ' num2str(pdx) ' ' ss ]  )
            mods(:,idx) = community_louvain2_squeezemore(clusterdat,prctile(triuvec(clusterdat,1),50),[],'potts',1e-6) ;
            disp(['resulted in n comms:' num2str(length(unique(mods(:,idx)))) ])
        end
          
        spkpttrnclust.(sdx{1}).mods{pdx} = mods ; 

    end
end

%% and get consensus

for sdx = subsets
        
    spkpttrnclust.(sdx{1}).cons = cell(3,1) ; 

    for pdx = 1:3
  
        disp([ 'consensus pdx' num2str(pdx)])

        agr =  agreement(spkpttrnclust.(sdx{1}).mods{pdx})./size(spkpttrnclust.(sdx{1}).mods{pdx},2) ; 
    
        spkpttrnclust.(sdx{1}).cons{pdx} = consensus_und2(agr,0.5,100) ; 

        disp(['resulted in n comms:' num2str(length(unique(spkpttrnclust.(sdx{1}).cons{pdx}))) ])

    end
end

% do a quick ordering of the consensus labels by size
for sdx = subsets
        
    spkpttrnclust.(sdx{1}).cons2 = cell(3,1)  ; 

    for pdx = 1:3
  
    
        ca = spkpttrnclust.(sdx{1}).cons{pdx} ; 

        uu = unique(ca) ; 
        for idx = 1:length(uu)
            ca(ca==uu(idx)) = idx ; 
        end

        tt = tabulate(ca) ; 
        %tt = tt(tt(:,2)>0,:) ;

        [~,szsortind] = sort(tt(:,2),'descend') ; 
        
        ca2 = zeros(length(ca),1) ; 
        for idx = 1:length(szsortind)
            ca2(ca==szsortind(idx)) = idx ; 
        end

        spkpttrnclust.(sdx{1}).cons2{pdx} = ca2 ; 

    end
end

%% now save all of this

filename = [ DD.PROC '/spk_rsspeaksclust_' OUTSTR '.mat' ] ; 
save(filename,'spkpttrnclust')

%% and analyze all the data we now have

