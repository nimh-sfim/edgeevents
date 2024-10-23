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

% fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

refparc = parc.ca(1:200) ; 
nedges = ( finfo.nnodes * (finfo.nnodes-1) ) / 2 ; 

%% load in motion and physio

loadConf = load_hcp_regressors('/Users/faskowitzji/joshstuff/data/hcp352_regressors/hcp352_regressors',...
      ['rfMRI_' imglob],'*rfMRI_REST1_RL_confounds.tsv') ; 

loadCC_WM = load_hcp_regressors('/Users/faskowitzji/joshstuff/data/hcp352_regressors/hcp352_regressors',...
      ['rfMRI_' imglob],'*rfMRI_REST1_RL_WM_acompcor.tsv') ; 

loadCC_CSF = load_hcp_regressors('/Users/faskowitzji/joshstuff/data/hcp352_regressors/hcp352_regressors',...
      ['rfMRI_' imglob],'*rfMRI_REST1_RL_CSF_acompcor.tsv') ; 

%% take the fft of the compcor components 

sdx = subsets ; 

fftWM = zeros(5,(finfo.ntp/2)+1) ; 
fftCSF = zeros(5,(finfo.ntp/2)+1) ; 

fftWM2 = zeros(5,129) ; 
fftCSF2 = zeros(5,129) ; 


for idx = 1:length(sublist.(sdx{1}))

    disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 

    % sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    % 
    % filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
    % readdat = load(filename,'spike_len_mat') ; 
    % 
    % ets = get_ets(datStr(sind).ts(:,1:finfo.nnodes)) ; 
    % 
    % dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
    % dd(isnan(dd)) = 0 ;
    % 
    % [~,o1] = count_spks(dd==1) ; 
    % [~,o2] = count_spks(dd==2) ; 
    % [~,o3] = count_spks(dd==3) ; 
    % 
    % etsrss_o = [ get_rss(o1) get_rss(o2) get_rss(o3) ] ; 
    % %etsrss_o = [ get_rss(dd==1) get_rss(dd==2) get_rss(dd==3) ] ; 


    %% physio signals 

    physigWM = table2array(loadCC_WM(datStr(sind).sub).values{1}) ;
    physigCSF = table2array(loadCC_CSF(datStr(sind).sub).values{1}) ;

    [~,ff ] = arrayfun(@(i_) quick_pspec(physigWM(:,i_),0.72),1:5,'UniformOutput',false) ;
    fftWM = cell2mat(ff)' + fftWM ; 

    [~,ff ] = arrayfun(@(i_) quick_pspec(physigCSF(:,i_),0.72),1:5,'UniformOutput',false) ;
    fftCSF = cell2mat(ff)' + fftCSF ; 


    [ff ] = arrayfun(@(i_) pwelch(physigWM(:,i_),20,0,[],1/0.72),1:5,'UniformOutput',false) ;
    fftWM2 = cell2mat(ff)' + fftWM2 ; 

    [ff ] = arrayfun(@(i_) pwelch(physigCSF(:,i_),20,0,[],1/0.72),1:5,'UniformOutput',false) ;
    fftCSF2 = cell2mat(ff)' + fftCSF2 ; 


end

fftWM = fftWM ./ length(sublist.(sdx{1})) ; 
fftCSF = fftCSF ./ length(sublist.(sdx{1})) ; 
fftWM2 = fftWM2 ./ length(sublist.(sdx{1})) ; 
fftCSF2 = fftCSF2 ./ length(sublist.(sdx{1})) ; 

% get the freq from pwelch
[~,f ] = pwelch(physigCSF(:,1),100,0,[],1/0.72) ; 

%% confirm that its 1 & 2nd comps that are respiration

cm = parula(5) ; 

TL = tiledlayout(1,2) ; 
TL.Title.String = 'pwelch aCompCor components'

nexttile

for idx = 1:5
plot(f,10*log10(fftWM2(idx,:)),'Color',cm(idx,:),'LineWidth',3)
hold on
end
hold off
title('aCompCor WM')

clim([0 1/0.72/2])

xlabel('freq (Hz)')
ylabel('PSD (dB/Hz)')

nexttile

for idx = 1:5
plot(f,10*log10(fftCSF2(idx,:)),'Color',cm(idx,:),'LineWidth',3)
hold on
end
hold off
title('aCompCor CSF')

xlabel('freq (Hz)')
ylabel('PSD (dB/Hz)')

clim([0 1/0.72/2])

lg = legend(cellstr([ repmat('comp. ',5,1) num2str((1:5)') ])) ; 
lg.Location = "northeastoutside"

set(gcf,'Position',[100 100 800 400])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/compcor_psd.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% now read in spike lengths to make a histogram of lengths

for sdx = subsets

    corr_w_phys.(sdx{1}).ccWMo = zeros(3,5,length(sublist.(sdx{1}))) ; 
    corr_w_phys.(sdx{1}).ccCSFo = zeros(3,5,length(sublist.(sdx{1}))) ; 

    corr_w_phys.(sdx{1}).xcov_long_WMo = zeros(85,length(sublist.(sdx{1}))) ; 
    corr_w_phys.(sdx{1}).xcov_long_CSFo = zeros(85,length(sublist.(sdx{1}))) ; 

    corr_w_phys.(sdx{1}).xcov_short_WMo = zeros(85,length(sublist.(sdx{1}))) ; 
    corr_w_phys.(sdx{1}).xcov_short_CSFo = zeros(85,length(sublist.(sdx{1}))) ; 


    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        ets = get_ets(datStr(sind).ts(:,1:finfo.nnodes)) ; 
            
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        [~,o1] = count_spks(dd==1) ; 
        [~,o2] = count_spks(dd==2) ; 
        [~,o3] = count_spks(dd==3) ; 

        etsrss_o = [ get_rss(o1) get_rss(o2) get_rss(o3) ] ; 
        %etsrss_o = [ get_rss(dd==1) get_rss(dd==2) get_rss(dd==3) ] ; 


        %% physio signals 

        physigWM = table2array(loadCC_WM(datStr(sind).sub).values{1}) ;
        physigCSF = table2array(loadCC_CSF(datStr(sind).sub).values{1}) ;

        %% corr

        corr_w_phys.(sdx{1}).ccWMo(:,:,idx) = (corr(etsrss_o,physigWM)) ; 
        corr_w_phys.(sdx{1}).ccCSFo(:,:,idx) = (corr(etsrss_o,physigCSF)) ; 

        %% xcov

        corr_w_phys.(sdx{1}).xcov_short_WMo(:,idx) = xcov(etsrss_o(:,1),physigWM(:,2),42,'normalized') ;
        corr_w_phys.(sdx{1}).xcov_short_CSFo(:,idx) = xcov(etsrss_o(:,1),physigCSF(:,2),42,'normalized') ;

        corr_w_phys.(sdx{1}).xcov_long_WMo(:,idx) = xcov(etsrss_o(:,3),physigWM(:,2),42,'normalized') ;
        corr_w_phys.(sdx{1}).xcov_long_CSFo(:,idx) = xcov(etsrss_o(:,3),physigCSF(:,2),42,'normalized') ;


    end
end



%% 

filename = [ DD.PROC '/spk_corrwphys_' OUTSTR '.mat' ] ; 
save(filename,'corr_w_phys','-v7.3')

%%

tiledlayout(2,2)

cc = inferno(10) ; 

plotnames = {'xcov_short_WMo' 'xcov_short_CSFo' 'xcov_long_WMo' 'xcov_long_CSFo' } ; 

longerplotnames1 = {'short' 'short' 'long' 'long'} ; 
longerplotnames2 = {'WM second comp.' 'CSF second comp.' 'WM second comp.' 'CSF second comp.'  } ; 

for idx = 1:4

    nexttile
    
    %plot(-42:1:42,corr_w_phys.subset1.xcovWM,'Color',[ cc(3,:) 0.2],'LineWidth',2)
    dat = corr_w_phys.subset1.(plotnames{idx}) ; 
    for ss = 1:size(dat,2)
        % plot(-42:1:42,dat(:,idx),'LineWidth',2)
        x = (-42:1:42) ; 
        y = dat(:,ss)' ; 
        z = zeros(size(dat(:,ss)))' ; 
        c = gradient(dat(:,ss))' ; 
        surface([x;x],[y;y],[z;z],[c;c],...
            'FaceColor','no','EdgeColor','interp','LineWidth',2,'EdgeAlpha',0.2)
        hold on
    end
    hold off
    % colormap(rdbu)
    colormap(flipud(inferno))
    xlim([-42 42])
    xline(0,'Color',[0.8 0.8 0.8 0.95],'LineWidth',4)
    ylim([-.5 .5])
    % and plot the mean
    
    hold on
    plot((-42:1:42),mean(dat,2),'Color',[0.5 0.5 0.5],'LineWidth',5)
    hold off
    
    xx = xticks ; 
    xticklabels(strtrim(cellstr(num2str(xx'.*finfo.TR))))
    xtickangle(0)

    clim([ -0.05 0.05 ])
    cb = colorbar() ; 
    cb.Label.String = 'derivative' ; 
    
    title([longerplotnames1{idx} ' vs. ' longerplotnames2{idx} ])

    if idx == 3 || idx == 4 
        xlabel('time (sec)')
    end

    if idx == 1 || idx == 3
        ylabel('correlation')
    end

end

set(gcf,'Position',[100 100 1000 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/xcorr_w_compcor.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

