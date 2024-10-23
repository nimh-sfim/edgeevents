%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%% make a tsnr map

subsets = {'subset1' 'subset2'} ; 

acfs = struct() ; 
for sdx = subsets

    acfs.(sdx{1}).map = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        tmp = zscore(datStr(sind).ts(:,1:finfo.nnodes)) ;

        acfs.(sdx{1}).map(:,idx) = get_acf_hwhm(tmp,0.72,20) ; 

    end

end

%%

for sdx = subsets

    acfs.(sdx{1}).mat = zeros(finfo.nnodes,finfo.nnodes) ; 

    for idx = 1:length(sublist.(sdx{1}))
    

        t = acfs.(sdx{1}).map(:,idx) ;
        acfs.(sdx{1}).mat = t*t' + acfs.(sdx{1}).mat ; 

    end

    acfs.(sdx{1}).mat = acfs.(sdx{1}).mat ./ length(sublist.(sdx{1})) ; 

end

%%

save('./data/interim/ts_autocorr.mat',"acfs")

%% load spk conn

filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
load(filename)

% nothing here
scatter_w_rho(tv(acfs.subset1.mat),tv(spike_conn.subset1.long))
