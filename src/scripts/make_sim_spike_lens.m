%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

maxspk = 100 ; 

high_bin = 4 ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

%% read in all the spikes

subsets = {'subset1' 'subset2'} ; 

% for sdx = subsets
% 
%     spike_lengths.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
% 
%     for idx = 1:length(sublist.(sdx{1}))
% 
%         disp(idx)
% 
%         sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
% 
%         filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
%         readdat = load(filename,'spike_len_cell') ; 
% 
%         spike_lengths.(sdx{1}){idx} = readdat.spike_len_cell  ; 
% 
%     end
% 
% end
% 
% for sdx = subsets
%     tmp = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.(sdx{1}){i_}')),1:length(spike_lengths.(sdx{1}))','UniformOutput',false)') ;
%     % get rid of 0 lengths
%     spike_lengths.all.(sdx{1}) = nonzeros(tmp) ; 
% end

%%

nperms = 1000 ;
simlens.keepcov = cell(nperms,1) ; 
simlens.nocov = cell(nperms,1) ; 
simlens.randcov = cell(nperms,1) ; 

simmat.keepcov = cell(nperms,1) ; 
simmat.nocov = cell(nperms,1) ; 
simmat.randcov = cell(nperms,1) ; 

simrssO.keepcov = cell(nperms,1) ; 
simrssO.nocov = cell(nperms,1) ; 
simrssO.randcov = cell(nperms,1) ; 

simrssA.keepcov = cell(nperms,1) ; 
simrssA.nocov = cell(nperms,1) ; 
simrssA.randcov = cell(nperms,1) ; 

%%

for idx = 1:nperms

    disp(['perm ' num2str(idx) ' of ' num2str(nperms)])

    % randomly pick cov to simulate
    randselect = randi(NSUBS); 

    ts = zscore(datStr(randselect).ts(:,1:finfo.nnodes)) ; 
    % tmpets = get_ets(ts) ; 
    % abv_thr = tmpets > SPK_THR ; 
    % [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ;
    % dd = discretize(spike_len_mat,lowmedhigh_edges) ;

    %% generate dat
    [simts_keepcov,simts_nocov] = simulate_BOLD_timecourse_func_v3(1200,0.72,0.72,cov(ts),mean((abs(fft(zscore(ts))).^2),2)) ; 

    % no cov

    tmpets = get_ets(simts_nocov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ;
    simlens.nocov{idx} = uint8(cell2mat(spike_len_cell')) ; 

    dd = discretize(spike_len_mat,lowmedhigh_edges) ;
    [c1,o1] = count_spks(dd==1) ; 
    [c2,o2] = count_spks(dd==2) ; 
    [c3,o3] = count_spks(dd==3) ; 

    simmat.nocov{idx} = [ c1 ; c2 ; c3 ] ; 
    simrssO.nocov{idx} = [get_rss(o1) get_rss(o2) get_rss(o3) ] ; 
    simrssA.nocov{idx} = [get_rss(dd==1) get_rss(dd==2) get_rss(dd==3) ] ;

    % keep cov

    tmpets = get_ets(simts_keepcov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ;
    simlens.keepcov{idx} = uint8(cell2mat(spike_len_cell')) ; % keep memory low

    dd = discretize(spike_len_mat,lowmedhigh_edges) ;
    [c1,o1] = count_spks(dd==1) ; 
    [c2,o2] = count_spks(dd==2) ; 
    [c3,o3] = count_spks(dd==3) ; 

    simmat.keepcov{idx} = [ c1 ; c2 ; c3 ] ; 
    simrssO.keepcov{idx} = [get_rss(o1) get_rss(o2) get_rss(o3) ] ; 
    simrssA.keepcov{idx} = [get_rss(dd==1) get_rss(dd==2) get_rss(dd==3) ] ;
    

    %% randcov 

    % [~,ev] = eig(cov(ts)) ; 
    % [simts_randcov] = simulate_BOLD_timecourse_func_v3(1200,0.72,0.72,...
    %     nearestSPD( gallery('randcorr',diag(ev)) ) ,mean((abs(fft(zscore(ts))).^2),2)) ; 

    [simts_randcov] = simulate_BOLD_timecourse_func_v3(1200,0.72,0.72,...
        nearestSPD( hqs(cov(ts)) ) ,mean((abs(fft(zscore(ts))).^2),2)) ; 

    tmpets = get_ets(simts_randcov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ; 

    simlens.randcov{idx} = uint8(cell2mat(spike_len_cell')) ; 

    dd = discretize(spike_len_mat,lowmedhigh_edges) ;
    [c1,o1] = count_spks(dd==1) ; 
    [c2,o2] = count_spks(dd==2) ; 
    [c3,o3] = count_spks(dd==3) ; 

    simmat.randcov{idx} = [ c1 ; c2 ; c3 ] ; 
    simrssO.randcov{idx} = [get_rss(o1) get_rss(o2) get_rss(o3) ] ; 
    simrssA.randcov{idx} = [get_rss(dd==1) get_rss(dd==2) get_rss(dd==3) ] ;
    

end

% save it

disp('SAVING!!')

filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spk.mat'] ; 
save(filename,'simlens','-v7.3')

filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
save(filename,'simmat','-v7.3')

filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spkrss.mat'] ; 
save(filename,'simrss*','-v7.3')


%%