%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2 ; 

for idx = 1:NSUBS

    disp(idx)

    filename = [DD.PROC '/' datStr(idx).sub '_' OUTSTR '_spike_len.mat'] ; 

    if isfile(filename)
        disp(['already finsiehd sub: ' datStr(idx).sub ])
        continue
    end

    tmpts = datStr(idx).ts ; 
    tmpets = get_ets(tmpts) ; 

    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ; 

    save([DD.PROC '/' datStr(idx).sub '_' OUTSTR '_spike_len.mat'],...
        'spike_len_cell','spike_len_mat','-v7.3')

end

%% now read in spike lengths to make a histogram of lengths

spike_lengths.subset1 = cell(length(sublist.subset1),1) ; 
spike_lengths.subset2 = cell(length(sublist.subset2),1) ; 

for idx = 1:length(sublist.subset1)

    disp(idx)

    sind = find(cellfun(@(x_)strcmp(x_,sublist.subset1(idx)),sublist.all)) ; 

    filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_spike_len.mat'] ; 
    readdat = load(filename,'spike_len_cell') ; 

    spike_lengths.subset1{idx} = readdat.spike_len_cell  ; 

end

all_spike_lengths = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.subset1{i_}')),1:length(spike_lengths.subset1)','UniformOutput',false)') ;

