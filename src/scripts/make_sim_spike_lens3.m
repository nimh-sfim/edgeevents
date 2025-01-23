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


%%

nperms = 1000 ;
simlens.keepcov = cell(nperms,1) ; 
simmat.keepcov = cell(nperms,1) ; 
simvar.keepcov = cell(nperms,1) ; 
simrssO.keepcov = cell(nperms,1) ; 
simrssA.keepcov = cell(nperms,1) ; 

simlens.nocov = cell(nperms,1) ; 
simmat.nocov = cell(nperms,1) ; 
simvar.nocov = cell(nperms,1) ; 
simrssO.nocov = cell(nperms,1) ; 
simrssA.nocov = cell(nperms,1) ; 

%%

for idx = 1:nperms

    disp(['perm ' num2str(idx) ' of ' num2str(nperms)])

    % randomly pick cov to simulate
    randselect = randi(NSUBS); 

    ts = zscore(datStr(randselect).ts(:,1:finfo.nnodes)) ; 

    %% generate dat
    % [simts_keepcov,~] = simulate_BOLD_timecourse_func_v3( ...
    %     finfo.ntp,finfo.TR,finfo.TR,cov(ts),mean((abs(fft(zscore(ts))).^2),2)) ; 
    [simts_nocov] = generate_phase_surrogates(ts,0,0) ; 

    tmpets = get_ets(simts_nocov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ;
    simlens.nocov{idx} = countcount(cell2mat(spike_len_cell'),1:100) ; % keep memory low

    simvar.nocov{idx} = cellfun(@(x_) mad(x_,0),spike_len_cell) ; 

    dd = discretize(spike_len_mat,lowmedhigh_edges) ;
    [c1,o1] = count_spks(dd==1) ; 
    [c2,o2] = count_spks(dd==2) ; 
    [c3,o3] = count_spks(dd==3) ; 
    [c4,o4] = count_spks(spike_len_mat >= 8) ; 
    [c5,o5] = count_spks(spike_len_mat >= 14) ; 

    simmat.nocov{idx} = [ c1 ; c2 ; c3 ; c4 ; c5 ] ; 
    simrssO.nocov{idx} = [sum(o1,2) sum(o2,2) sum(o3,2) sum(o4,2) sum(o5,2) ] ; 
    simrssA.nocov{idx} = [sum(dd==1,2) sum(dd==2,2) sum(dd==3,2) ] ;
    
    % keep cov

    [simts_keepcov] = generate_phase_surrogates(ts,1,0) ; 

    tmpets = get_ets(simts_keepcov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ;
    simlens.keepcov{idx} = countcount(cell2mat(spike_len_cell'),1:100) ; % keep memory low ; % keep memory low

    simvar.keepcov{idx} = cellfun(@(x_) mad(x_,0),spike_len_cell) ; 

    dd = discretize(spike_len_mat,lowmedhigh_edges) ;
    [c1,o1] = count_spks(dd==1) ; 
    [c2,o2] = count_spks(dd==2) ; 
    [c3,o3] = count_spks(dd==3) ; 
    [c4,o4] = count_spks(spike_len_mat >= 8) ; 
    [c5,o5] = count_spks(spike_len_mat >= 14) ; 

    simmat.keepcov{idx} = [ c1 ; c2 ; c3 ; c4 ; c5 ] ; 
    simrssO.keepcov{idx} = [sum(o1,2) sum(o2,2) sum(o3,2) sum(o4,2) sum(o5,2)] ; 
    simrssA.keepcov{idx} = [sum(dd==1,2) sum(dd==2,2) sum(dd==3,2) ] ;
    
end

% save it

disp('SAVING!!')

filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spk.mat'] ; 
save(filename,'simlens','-v7.3')

filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
save(filename,'simmat','simvar','-v7.3')

filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spkrss.mat'] ; 
save(filename,'simrss*','-v7.3')


%%

mats = zeros(200,200,3) ; 

for idx = 1:1000
    disp(idx)

    for jdx = 1:3

        mats(:,:,jdx) =  mats(:,:,jdx) + mksq(mm.simmat.nocov{idx}(jdx,:)) ; 
    end
end

% load spk conn
filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
load(filename,'spike_conn')
 