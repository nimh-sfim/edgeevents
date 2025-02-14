%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 
freqStr = struct(); 

for sdx = subsets

    tmp1 = zeros(557,length(sublist.(sdx{1}))) ; 
    tmp2 = zeros(557,length(sublist.(sdx{1}))) ; 

    parfor idx = 1:length(sublist.(sdx{1}))
    
        disp([num2str(idx) ' - ' sdx{1}])
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        ts = datStr(sind).ts(:,1:finfo.nnodes) ; 
        ets = get_ets(ts) ; 

        [pp1] = pspectrum(ts,1/0.72,'FrequencyResolution',0.01) ; 
        [pp2] = pspectrum(single(ets>SPK_THR),1/0.72,'FrequencyResolution',0.01) ;
        
        tmp1(:,idx) = mean(pp1,2,'omitmissing') ; 
        tmp2(:,idx) = mean(pp2,2,'omitmissing') ; 

    end


    freqStr.(sdx{1}).ts = tmp1 ; 
    freqStr.(sdx{1}).spk = tmp2 ; 

end

ts = datStr(42).ts(:,1:finfo.nnodes) ; 
[~,sampFreqs] = pspectrum(ts,1/0.72,'FrequencyResolution',0.01) ; 

%%

filename = [ DD.PROC '/freq_info_' OUTSTR '.mat' ] ; 
save(filename,'freqStr','sampFreqs', '-v7.3')

%%  make a plot!

tiledlayout(2,2)

nexttile()
plot_manylines(sampFreqs,log10(freqStr.subset1.ts)*100,'Color',[0.2 0.2 0.2 0.2])
xlim([min(sampFreqs) max(sampFreqs)])

xlabel('frequency')
ylabel('power (dB)')

title('time series')

nexttile()
plot_manylines(sampFreqs,log10(freqStr.subset1.spk)*100,'Color',[0.2 0.2 0.2 0.2])
xlim([min(sampFreqs) max(sampFreqs)])

xlabel('frequency')
ylabel('power (dB)')

title('edge events')

nexttile()
plot_manylines(sampFreqs,log10(freqStr.subset2.ts)*100,'Color',[0.2 0.2 0.2 0.2])
xlim([min(sampFreqs) max(sampFreqs)])

xlabel('frequency')
ylabel('power (dB)')


nexttile()
plot_manylines(sampFreqs,log10(freqStr.subset2.spk)*100,'Color',[0.2 0.2 0.2 0.2])
xlim([min(sampFreqs) max(sampFreqs)])

xlabel('frequency')
ylabel('power (dB)')


set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/across_subject_pspec.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)