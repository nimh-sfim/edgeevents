%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

%% load the spike conn

filename = [ DD.PROC '/spk_conn_indiv_' OUTSTR '.mat' ] ; 
l1 = load(filename) ; 

filename = [ DD.PROC '/spk_conn_indiv_' OUTSTR '_2.mat' ] ; 
l2 = load(filename) ; 

filename = [ DD.PROC '/spk_conn_ex_indiv_' OUTSTR '.mat' ] ; 
l1ex = load(filename) ; 

filename = [ DD.PROC '/spk_conn_ex_indiv_' OUTSTR '_2.mat' ] ; 
l2ex = load(filename) ; 

%% identifiability w/ spk counts

ident = struct() ; 
identpos = struct() ; 

for sdx = l1.spklen_names

    ident.(sdx{1}) = zeros(176) ; 
    
    for idx = 1:length(sublist.subset1)
        disp(idx)
        c1 = (l1.spike_conn_indiv.subset1.(sdx{1})(:,:,idx)) ; 
        c1 = c1(1:finfo.nnodes,1:finfo.nnodes) ; 
    
        for jdx = 1:length(sublist.subset1)
            c2 = (l2.spike_conn_indiv.subset1.(sdx{1})(:,:,jdx)) ; 
            c2 = c2(1:finfo.nnodes,1:finfo.nnodes) ; 
           
            ident.(sdx{1})(idx,jdx) = IPN_ccc([ tv(c1) tv(c2) ]) ; 
            % ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
            % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
            % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
            % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
        end
    end

    [~,ss1] = sort(ident.(sdx{1}),2,'descend') ; 
    identpos.(sdx{1}) = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 
    
end

for sdx = {'longest'}

    ident.(sdx{1}) = zeros(176) ; 
    
    for idx = 1:length(sublist.subset1)
        disp(idx)
        c1 = mksq(l1ex.spike_conn_ex_indiv.subset1.(sdx{1})(:,idx)) ; 
        c1 = c1(1:finfo.nnodes,1:finfo.nnodes) ; 
    
        for jdx = 1:length(sublist.subset1)
            c2 = mksq(l2ex.spike_conn_ex_indiv.subset1.(sdx{1})(:,jdx)) ; 
            c2 = c2(1:finfo.nnodes,1:finfo.nnodes) ; 
           
            ident.(sdx{1})(idx,jdx) = IPN_ccc([ tv(c1) tv(c2) ]) ; 
            % ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
            % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
            % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
            % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
        end
    end

    [~,ss1] = sort(ident.(sdx{1}),2,'descend') ; 
    identpos.(sdx{1}) = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 
    
end


% and finally do with staight up FC
filename = [ DD.PROC '/fc_conn_indiv_' OUTSTR '.mat' ] ; 
f1 = load(filename) ; 

filename = [ DD.PROC '/fc_conn_indiv_' OUTSTR '_2.mat' ] ; 
f2 = load(filename) ; 

for sdx = {'fc'}

    ident.(sdx{1}) = zeros(176) ; 

    for idx = 1:length(sublist.subset1)
        disp(idx)
        c1 = (f1.fc_conn_indiv.subset1(:,:,idx)) ; 
        c1 = c1(1:finfo.nnodes,1:finfo.nnodes) ; 

        for jdx = 1:length(sublist.subset1)
            c2 = (f2.fc_conn_indiv.subset1(:,:,jdx)) ; 
            c2 = c2(1:finfo.nnodes,1:finfo.nnodes) ; 

            ident.(sdx{1})(idx,jdx) = corr(tv(c1),tv(c2),'type','p') ; 
            % ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
            % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
            % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
            % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
        end
    end

    [~,ss1] = sort(ident.(sdx{1}),2,'descend') ; 
    identpos.(sdx{1}) = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 

end

%%

% scannames = {'REST1_LR' 'REST1_RL'} ; 
% 
% fmridat = struct() ; 
% for idx = scannames
%     fmridat.(idx{1}) = load_hcp_alldata(...
%         [ DD.INTERM '/hcp352_nusregts_FIX2phys_schaefer200/' ],...
%         'schaefer200-yeo17',...
%         sublist.subset1, ...
%         ['*' idx{1} '*']) ; 
% end


%%

filename = [ DD.PROC '/spk_longestedge_sigmask_' OUTSTR '.mat' ] ; 
load(filename)

dynident = struct() ; 
dynident.longest = zeros(176) ; 
% dynidentpos = struct() ; 

sigmask = logical(sig.B.more.longest) ; 

%%
% 
% dynmats.REST1_RL.mats = nan(finfo.ntp*( finfo.ntp-1)/2,176) ; 
% dynmats.REST1_LR.mats = nan(finfo.ntp*( finfo.ntp-1)/2,176) ; 
% 
% for idx = 1:176
%     s1 = fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes) ; 
%     e1 = get_ets(s1) ; 
%     d1 = corr(e1(:,tv(sigmask))') ;
%     dynmats.REST1_RL.mats = tv() ; 
% end

%%

% for idx = 13:length(sublist.subset1)
%     disp(idx)
%     s1 = fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes) ; 
%     e1 = get_ets(s1) ; 
%     d1 = corr(e1(:,tv(sigmask))') ; 
% 
%     for jdx = 1:length(sublist.subset1)
%         s2 = fmridat.REST1_LR(jdx).ts(:,1:finfo.nnodes) ; 
%         e2 = get_ets(s2) ; 
%         d2 = corr(e2(:,tv(sigmask))') ; 
% 
%         dynident.longest(idx,jdx) =  adj_spect_dist(d1,d2,10) ; 
%         % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
%         % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
%         % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
%     end
% end

% [~,mi1] = max(ident.long,[],2) ; 
% [~,mi2] = max(ident.notlong,[],2) ; 
% 

% %% now read in spike lengths 
% 
% subsets = {'subset1'} ; 
% 
% for sdx = subsets
% 
%     spike_mats.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
% 
%     for idx = 1:length(sublist.(sdx{1}))
% 
%         disp(idx)
% 
%         sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
% 
%         filename = [DD.PROC '/' 'REST1_RL' '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
%         readdat = load(filename,'spike_len_mat') ; 
% 
%         spike_mats.REST1_RL.(sdx{1}){idx} = readdat.spike_len_mat  ; 
% 
%         filename = [DD.PROC '/' 'REST1_LR' '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
%         readdat = load(filename,'spike_len_mat') ; 
% 
%         spike_mats.REST1_LR.(sdx{1}){idx} = readdat.spike_len_mat  ; 
% 
%     end
% 
% end

%%

% maxspk = 100 ; 
% 
% high_bin = 4 ; 
% lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 
% 
% ident2 = struct() ; 
% identpos2 = struct() ; 
% 
% % lets first male all the disc data
% 
% for idx = 1:length(sublist.subset1)
%     disp(idx)
% 
%     dd = discretize(spike_mats.REST1_RL.subset1{idx} ,lowmedhigh_edges) ; 
%     dd(isnan(dd)) = 0 ; 
%     spike_mats.REST1_RL.subset1{idx} = dd ; 
% 
%     dd = discretize(spike_mats.REST1_LR.subset1{idx} ,lowmedhigh_edges) ; 
%     dd(isnan(dd)) = 0 ; 
%     spike_mats.REST1_LR.subset1{idx} = dd ; 
% 
% end
% 
% spkts = struct() ; 
% for sdx = 1:length(l1.spklen_names)
% 
%     spkts.REST1_RL.subset1.(l1.spklen_names{sdx}) = cell(176,1) ; 
%     spkts.REST1_LR.subset1.(l1.spklen_names{sdx}) = cell(176,1) ; 
% 
%     for idx = 1:length(sublist.subset1)
%         disp(idx)
% 
%         spkts.REST1_RL.subset1.(l1.spklen_names{sdx}){idx} = ets_2_ts(spike_mats.REST1_RL.subset1{idx}==sdx) ; 
%         spkts.REST1_LR.subset1.(l1.spklen_names{sdx}){idx} = ets_2_ts(spike_mats.REST1_LR.subset1{idx}==sdx) ; 
% 
% 
%     end
% end
% 
% %%
% 
% ident2 = struct() ; 
% identpos2 = struct() ; 
% 
% for sss = 1:length(l1.spklen_names)
% 
%     sdx = l1.spklen_names(sss) ; 
% 
%     ident2.(sdx{1}) = zeros(176) ; 
% 
%     for idx = 1:length(sublist.subset1)
%         disp(idx)
%         e1 = spkts.REST1_RL.subset1.(sdx{1}){idx} ; 
%         c1 = corr(e1) ; 
%         c1 = c1(1:finfo.nnodes,1:finfo.nnodes) ; 
% 
%         for jdx = 1:length(sublist.subset1)
%             %disp(jdx)
%             e2 = spkts.REST1_LR.subset1.(sdx{1}){jdx} ; 
%             c2 = corr(e2) ; 
%             c2 = c2(1:finfo.nnodes,1:finfo.nnodes) ; 
% 
%             ident2.(sdx{1})(idx,jdx) = corr(tv(c1),tv(c2),'type','p') ; 
%             % ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
%             % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
%             % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
%             % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
%         end
%     end
% 
%     [~,ss1] = sort(ident2.(sdx{1}),2,'descend') ; 
%     identpos2.(sdx{1}) = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 
% 
% end
