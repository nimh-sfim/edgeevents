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

%% now read in spike lengths to make a histogram of lengths

for sdx = subsets(1)

    corr_w_mot.(sdx{1}).ccMoto = zeros(3,length(sublist.(sdx{1}))) ; 
    corr_w_mot.(sdx{1}).xcov_inter_Moto = zeros(85,length(sublist.(sdx{1}))) ; 
    corr_w_mot.(sdx{1}).xcov_long_Moto = zeros(85,length(sublist.(sdx{1}))) ; 
    corr_w_mot.(sdx{1}).xcov_short_Moto = zeros(85,length(sublist.(sdx{1}))) ; 

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

        %% physio signals 

        motsig = loadConf(datStr(sind).sub).values{1}.framewise_displacement ;  

        %% corr

        corr_w_mot.(sdx{1}).ccMoto(:,idx) = (corr(etsrss_o,motsig)) ; 

        %% xcov

        corr_w_mot.(sdx{1}).xcov_short_Moto(:,idx) = xcov(etsrss_o(:,1),motsig,42,'normalized') ;
        corr_w_mot.(sdx{1}).xcov_inter_Moto(:,idx) = xcov(etsrss_o(:,2),motsig,42,'normalized') ;
        corr_w_mot.(sdx{1}).xcov_long_Moto(:,idx) = xcov(etsrss_o(:,3),motsig,42,'normalized') ;


    end
end


%% 

filename = [ DD.PROC '/spk_corrwmot_' OUTSTR '.mat' ] ; 
save(filename,'corr_w_mot','-v7.3')

%%

tiledlayout(1,3)

cc = inferno(10) ; 

plotnames = {'xcov_short_Moto' 'xcov_inter_Moto' 'xcov_long_Moto' } ; 

longerplotnames1 = {'short' 'inter' 'long' } ; 
longerplotnames2 = {'mot' 'mot' 'mot'  } ; 

for idx = 1:3

    nexttile
    
    %plot(-42:1:42,corr_w_phys.subset1.xcovWM,'Color',[ cc(3,:) 0.2],'LineWidth',2)
    dat = corr_w_mot.subset1.(plotnames{idx}) ; 
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

    if idx == 2
        xlabel('time (sec)')
    end

    if idx == 1 
        ylabel('correlation')
    end

end

set(gcf,'Position',[100 100 800 200])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/xcorr_w_mot.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

tiledlayout(1,3,'TileIndexing','rowmajor')

nexttile()
histogram(squeeze(corr_w_mot.(sdx{1}).ccMoto(1,:)))
title('short vs mot.')

ylabel('count')

nexttile()
histogram(squeeze(corr_w_mot.(sdx{1}).ccMoto(2,:)))
title('inter vs mot.')

nexttile()
histogram(squeeze(corr_w_mot.(sdx{1}).ccMoto(3,:)))
title('long vs mot.')


set(gcf,'Position',[100 100 800 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/corr_w_mot.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)




