%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

%% read in all the spikes

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_lengths.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
    
        spike_lengths.(sdx{1}){idx} = readdat.spike_len_cell  ; 
    
    end

end

for sdx = subsets
    tmp = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.(sdx{1}){i_}')),1:length(spike_lengths.(sdx{1}))','UniformOutput',false)') ;
    % get rid of 0 lengths
    spike_lengths.all.(sdx{1}) = nonzeros(tmp) ; 
end

%%

nperms = 500 ;
simlens_keepcov = cell(nperms,1) ; 
simlens_nocov = cell(nperms,1) ; 
simlens_randcov = cell(nperms,1) ; 

for idx = 1:nperms

    disp(['perm ' num2str(idx) ' of ' num2str(nperms)])

    % randomly pick cov to simulate
    randselect = randi(NSUBS); 

    ts = zscore(datStr(randselect).ts) ; 

    [simts_keepcov,simts_nocov] = simulate_BOLD_timecourse_func_v3(1200,0.72,0.72,cov(ts),mean((abs(fft(zscore(ts))).^2),2)) ; 

    % keep cov

    tmpets = get_ets(simts_keepcov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [~,spike_len_cell] = spk_lenmat(abv_thr) ; 

    simlens_keepcov{idx} = uint8(cell2mat(spike_len_cell')) ; % keep memory low

    % no cov

    tmpets = get_ets(simts_nocov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [~,spike_len_cell] = spk_lenmat(abv_thr) ; 

    simlens_nocov{idx} = uint8(cell2mat(spike_len_cell')) ; 

    % randcov 

    [simts_randcov] = simulate_BOLD_timecourse_func_v3(1200,0.72,0.72,...
        nearestSPD( hqs(cov(ts)) ) ,mean((abs(fft(zscore(ts))).^2),2)) ; 

    tmpets = get_ets(simts_randcov) ; 
    abv_thr = tmpets > SPK_THR ; 
    [~,spike_len_cell] = spk_lenmat(abv_thr) ; 

    simlens_randcov{idx} = uint8(cell2mat(spike_len_cell')) ; 

end

%% save it

filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spk.mat'] ; 
save(filename,'simlens_*','-v7.3')