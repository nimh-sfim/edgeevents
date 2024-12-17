%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_2.m') 

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_mats.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        spike_mats.(sdx{1}){idx} = readdat.spike_len_mat  ; 
    
    end

end

%% 

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

spk_amplencorr = struct() ; 

for sdx = subsets

    spk_amplencorr.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        % tmp_spkmat = spike_mats.(sdx{1}){idx}(:,cortmask) > 0 ;  
        tmp_ets = get_ets(datStr(sind).ts(:,1:finfo.nnodes)) ; 
        tmp_spk = tmp_ets>SPK_THR ; 

        % get peak amp per spike len
        [~,~,~,dd] = arrayfun(@(i_) get_contact_times(tmp_spk(:,i_)), 1:size(tmp_ets,2), 'UniformOutput', false) ; 

        % e1 = nan(size(tmp_ets,2),1) ; 
        % e2 = nan(size(tmp_ets,2),1) ; 
        cc = nan(size(tmp_ets,2),1) ; 
        for jdx = 1:size(tmp_ets,2)
            
            if isnan(dd{jdx})
                continue
            end

            e1  = arrayfun(@(i_) max(tmp_ets(dd{jdx}(i_,1):dd{jdx}(i_,2),jdx)),1:size(dd{jdx},1)) ; 
            e2  = arrayfun(@(i_) length(tmp_ets(dd{jdx}(i_,1):dd{jdx}(i_,2),jdx)),1:size(dd{jdx},1)) ; 

            cc(jdx) = corr(e1(:),e2(:),'type','s') ;  
        end
        cc(isnan(cc)) = 0 ; 

        spk_amplencorr.(sdx{1}) = spk_amplencorr.(sdx{1}) + mksq(cc) ; 

    end

    spk_amplencorr.(sdx{1}) = spk_amplencorr.(sdx{1}) ./ length(sublist.(sdx{1})) ; 
        
end

%% save this data

filename = [ DD.PROC '/spk_amplencorr_' OUTSTR '.mat' ] ; 
save(filename,'spk_amplencorr','-v7.3')
