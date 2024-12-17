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

%% that works

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
           
            ident.(sdx{1})(idx,jdx) = corr(tv(c1),tv(c2),'type','p') ; 
            % ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
            % ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
            % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
            % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
        end
    end
    % 
    % [~,mi1] = max(ident.long,[],2) ; 
    % [~,mi2] = max(ident.notlong,[],2) ; 
    
    [~,ss1] = sort(ident.(sdx{1}),2,'descend') ; 
    identpos.(sdx{1}) = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 
    
    % [~,ss2] = sort(ident.notlong,2,'descend') ; 
    % identpos2 = arrayfun(@(i_) find(ss2(i_,:)==i_),1:size(ss1,1)) ; 

end

%%

scannames = {'REST1_LR' 'REST1_RL'} ; 

fmridat = struct() ; 
for idx = scannames
    fmridat.(idx{1}) = load_hcp_alldata(...
        [ DD.INTERM '/hcp352_nusregts_FIX2phys_schaefer200/' ],...
        'schaefer200-yeo17',...
        sublist.subset1, ...
        ['*' idx{1} '*']) ; 
end


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

