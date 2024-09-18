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

fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

%% now read in spike lengths to make a histogram of lengths

for sdx = subsets

    spike_sys_trans.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ; 

        systrnStr = struct() ; 

        for jdx = 1:3

            % transition probability

            s1 = get_ets_states(single(dd==jdx),parc.ca(1:200),'assortweak') ; 

            [~,state_vec] = max(s1) ; 
            state_vec(sum(s1)==0) = 0 ; % correct for timepoints with not sys

            if length(unique(state_vec)) == 18
                tmp = make_trans_prob(state_vec,'uniq') ;  
            else
                warning('not all states detected')
                tmp = nan(18) ;
                inds = ismember(unique(state_vec),0:17) ; 
                tmp(inds,inds) = make_trans_prob(state_vec,'uniq') ; 
            end

            systrnStr.trans.(spklen_names{jdx}) = tmp ;

            % transition density
            tmp = cell2mat(arrayfun(@(i_) ...
                fullunroll(make_trans_dens(s1>0,i_,0,@mean)) , 2:21 , 'UniformOutput' , false)) ; 
            systrnStr.dens.(spklen_names{jdx}) = tmp ; 

            % transition onset density
            [~,s2] = count_spks((s1>0)') ; 
            tmp = cell2mat(arrayfun(@(i_) ...
                fullunroll(make_trans_dens(s2',i_,0,@any)) , 2:21 , 'UniformOutput' , false)) ; 
            systrnStr.tonst.(spklen_names{jdx}) = tmp ;      

        end

        spike_sys_trans.(sdx{1}){idx} = systrnStr  ; 
    
    end

end

%%

filename = [ DD.PROC '/spk_systrans_' OUTSTR '.mat' ] ; 
save(filename,'spike_sys_trans','-v7.3')

%% let's get some average mats

avg_sys_trans = struct() ; 
trans_names = {'trans' 'dens' 'tonst'} ; 

t_win = 15 ; 

for sdx = subsets
    
    for jdx = 1:3
        avg_sys_trans.trans.(sdx{1}).(spklen_names{jdx}) = zeros(17) ; 
        avg_sys_trans.dens.(sdx{1}).(spklen_names{jdx}) = zeros(17) ; 
        avg_sys_trans.tonst.(sdx{1}).(spklen_names{jdx}) = zeros(17) ; 

    end

    for idx = 1:length(sublist.(sdx{1}))

            for jdx = 1:3

                %% prob trans
        
                tmp1 = spike_sys_trans.(sdx{1}){idx}.trans.(spklen_names{jdx}) ; 
                tmp2 = tmp1(2:end,2:end) ; 
                tmp3 = normalize(tmp2,2,'norm',1) ; 
                tmp3(isnan(tmp3)) = 0 ; 

                avg_sys_trans.trans.(sdx{1}).(spklen_names{jdx}) = ...
                    tmp3 + ...
                    avg_sys_trans.trans.(sdx{1}).(spklen_names{jdx}) ; 

                %% dens

                tmp1 = spike_sys_trans.(sdx{1}){idx}.dens.(spklen_names{jdx}) ; 
                tmp2 = reshape(tmp1(:,t_win),17,17) ; 
                tmp3 = normalize(tmp2,2,'norm',1) ; % makes rows add to 1
                % tmp3 = tmp2 ./ diag(tmp2) ; % makes infs, not good idea
                tmp3(isnan(tmp3)) = 0 ; 

                avg_sys_trans.dens.(sdx{1}).(spklen_names{jdx}) = ...
                    tmp3 + ...
                    avg_sys_trans.dens.(sdx{1}).(spklen_names{jdx}) ; 
                
                %% onset

                tmp1 = spike_sys_trans.(sdx{1}){idx}.tonst.(spklen_names{jdx}) ; 
                tmp2 = reshape(tmp1(:,t_win),17,17) ; 
                tmp3 = normalize(tmp2,2,'norm',1) ; % makes rows add to 1
                tmp3(isnan(tmp3)) = 0 ; 

                avg_sys_trans.tonst.(sdx{1}).(spklen_names{jdx}) = ...
                    tmp3 + ...
                    avg_sys_trans.tonst.(sdx{1}).(spklen_names{jdx}) ; 

            end
    end

    for jdx = 1:3

        for tdx = 1:3

            avg_sys_trans.(trans_names{tdx}).(sdx{1}).(spklen_names{jdx}) = ...
                avg_sys_trans.(trans_names{tdx}).(sdx{1}).(spklen_names{jdx}) ./ length(sublist.(sdx{1})) ; 
        
        end

     end 
end

%% lets organize the systems based on meanfc pc?

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

% get mean values
[~,g1sort] = sort(arrayfun(@(i_) mean(grad_cifti(1,parc.ca(1:200)==i_)), 1:17)) ; 
[~,g2sort] = sort(arrayfun(@(i_) mean(grad_cifti(2,parc.ca(1:200)==i_)), 1:17)) ; 

% example viz
% imsc_grid_comm(avg_sys_trans.tonst.subset1.long(g1sort,g1sort),1:17,[],[],[],parc.names(g1sort)) ; clim([0 0.1])

