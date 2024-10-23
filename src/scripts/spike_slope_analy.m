%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

high_bin = 4 ; 
maxspk = 1200  ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

spklen_names = {'short' 'inter' 'long'} ; 

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

subsets = {'subset1' 'subset2'} ; 
nedges = ( finfo.nnodes * (finfo.nnodes-1) ) / 2 ; 

edgetrim = 42 ;
timevec = 1:finfo.ntp ;  
wanttime = timevec((edgetrim/2):end-(edgetrim/2)-1) ;

%% get the empirical slopes and the kendall

spike_slope = struct() ; 

for sdx = subsets

    % spike_slope.(sdx{1}).slope = zeros(3,length(sublist.(sdx{1}))) ; 
    spike_slope.(sdx{1}).slope_o = zeros(3,length(sublist.(sdx{1}))) ; 

    % spike_slope.(sdx{1}).rss = zeros(length(sublist.(sdx{1})),finfo.ntp-edgetrim,3) ; 
    spike_slope.(sdx{1}).rss_o = zeros(length(sublist.(sdx{1})),finfo.ntp-edgetrim,3) ; 

    % spike_slope.(sdx{1}).k = zeros(length(sublist.(sdx{1})),finfo.ntp-edgetrim,3) ; 
    spike_slope.(sdx{1}).k_o = zeros(3,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        %ets = get_ets(datStr(sind).ts(:,1:finfo.nnodes)) ; 
            
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        % origfc = corr(datStr(sind).ts(:,1:200)) ; 

        %% lets get the 3 rss 

        [~,o1] = count_spks(dd==1) ; 
        [~,o2] = count_spks(dd==2) ; 
        [~,o3] = count_spks(dd==3) ; 

        %  tmprss = [get_rss(dd==1) get_rss(dd==2) get_rss(dd==3)] ; 
        % tmprss_o = [get_rss(o1) get_rss(o2) get_rss(o3)] ; 
        % tmprss_o = tmprss_o(wanttime,:) ; 
        % change to sum cause it simple 
        tmprss_o = [sum(o1,2) sum(o2,2) sum(o3,2)] ; 
        tmprss_o = tmprss_o(wanttime,:) ; 

        % linear slope

        for jdx = 1:3
            % rr = regress(tmprss(:,jdx),[ ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)' ]) ;
            % spike_slope.(sdx{1}).slope(jdx,idx) = rr(2) ; 
                
            [rr,~] = regress(tmprss_o(:,jdx),[ ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)' ]) ;
            spike_slope.(sdx{1}).slope_o(jdx,idx) = rr(2) ; 
        end

        % rsss

        % spike_slope.(sdx{1}).rss(idx,:,1) = tmprss(:,1) ; 
        % spike_slope.(sdx{1}).rss(idx,:,2) = tmprss(:,2) ; 
        % spike_slope.(sdx{1}).rss(idx,:,3)  = tmprss(:,3) ; 

        spike_slope.(sdx{1}).rss_o(idx,:,1) = tmprss_o(:,1) ; 
        spike_slope.(sdx{1}).rss_o(idx,:,2) = tmprss_o(:,2) ; 
        spike_slope.(sdx{1}).rss_o(idx,:,3)  = tmprss_o(:,3) ; 

        % kendall

        for jdx = 1:3

            rr = Modified_MannKendall_test(...
                (1:finfo.ntp-edgetrim)',tmprss_o(:,jdx),0.05,0.05) ;
            spike_slope.(sdx{1}).k_o(jdx,idx) = rr ; 
        end
   

    end
end

%% save it

filename = [ DD.PROC '/spk_slope_' OUTSTR '.mat' ] ; 
save(filename,'spike_slope','-v7.3')

%% 

tiledlayout(3,1)

for idx = 1:3
    nexttile()

    [~,sind] = sort(spike_slope.(sdx{1}).slope_o(idx,:)) ; 

    tmp=spike_slope.(sdx{1}).rss_o(sind,:,idx) ; 

    imagesc(tmp,[0 100])
    
    %xlim([min(wanttime) max(wanttime)])
    yticks('')

    xx = xticklabels ; 
    xticklabels(num2str(round((str2double(xx)+(edgetrim/2)).*0.72,0)))

    if idx == 2
        ylabel('subjects')
        cb = colorbar ; 
        cb.Label.String = 'event onset count' ; 
    end

    if idx == 3
        xlabel('time (sec)')
    end

end

colormap(bone(100))

set(gcf,'Position',[100 100 800 400])
set(gcf,'Color','w')

%%


out_figdir = [ './reports/figures/figE/' ]
mkdir(out_figdir)
filename = [out_figdir '/slope_subjects_stacked.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)



%% do some stat testing

% load up the nulls
% filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spkrss.mat'] ; 
% load(filename)

%% bootstrap the mean slope, also bootstrap the regression coefs

nboot = 5000 ; 

for sdx = subsets

    spike_slope.(sdx{1}).boot_slopes = zeros(3,nboot) ; 

    subsz = length(sublist.(sdx{1})) ; 

    rng(42)
    for idx = 1:nboot

        disp(['boot: ' num2str(idx)])

        % select a random subset with replacement
        bootinds = randsample(1:subsz,subsz,true) ; 
 
        spike_slope.(sdx{1}).boot_slopes(:,idx) = ...
            mean(spike_slope.(sdx{1}).slope_o(:,bootinds),2) ; 
    end

end

%% bootstrap the regession coefs

aa = arrayfun(@(i_) regress((mean(spike_slope.(sdx{1}).rss_o(:,:,i_)))', [ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)']),1:3,'UniformOutput',false) ; 

for sdx = subsets

    spike_slope.(sdx{1}).boot_coefs = cell(nboot,1) ; 

    subsz = length(sublist.(sdx{1})) ; 

    rng(42)
    for idx = 1:nboot

        disp(['boot: ' num2str(idx)])

        % select a random subset with replacement
        bootinds = randsample(1:subsz,subsz,true) ; 
 
        rr = arrayfun(@(i_) ...
            regress((mean(spike_slope.(sdx{1}).rss_o(bootinds,:,i_)))', ...
            [ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)']), ...
            1:3,'UniformOutput',false) ; 

        spike_slope.(sdx{1}).boot_coefs{idx} = rr ; 
     
    end

end

%% make the plot with boot conf interval

cm = plasma(3) ; 

tiledlayout(3,2,"TileIndexing","columnmajor")

aa = arrayfun(@(i_) regress((mean(spike_slope.(sdx{1}).rss_o(:,:,i_)))', [ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)']),1:3,'UniformOutput',false) ; 

sdx = { 'subset1' } ; 

for idx = 1:3

    nexttile()

    xvals = (1:finfo.ntp-edgetrim) ; 

    % the data
    plot(xvals,(mean(spike_slope.(sdx{1}).rss_o(:,:,idx))),'Color',[0.2 0.2 0.2])
    hold on

    bootlines = zeros(finfo.ntp-edgetrim,nboot) ; 
    xx = [ones(finfo.ntp-edgetrim,1) (1:finfo.ntp-edgetrim)'] ; 
    for bdx = 1:nboot
        bootlines(:,bdx) = xx*spike_slope.(sdx{1}).boot_coefs{bdx}{idx} ; 
    end
    bootint = prctile(bootlines,[1 99],2) ; 
    fill([xvals fliplr(xvals)],[bootint(:,1)' flipud(bootint(:,2))'], ...
        [ 0.75 0.75 0.75 ],'facealpha',.5,'edgealpha',0,'FaceColor',cm(idx,:));
    xlim([xvals(1) xvals(end)])

    % the full fit
    rr = regress(mean(spike_slope.(sdx{1}).rss_o(:,:,idx))',xx) ; 
    plot(xvals,xx*rr,'Color',cm(idx,:),'LineWidth',2)
    hold off

    if idx == 3
        xl = xticks ; 
        xticklabels(round((xl+edgetrim).*finfo.TR))
        xlabel('time (sec)')
    else
        xticklabels('')
    end
    
    if idx == 2
        ylabel('num. event onsets')
    end
    

end

% the histograms of the slopes 
nexttile([3 1])

meanslopes = mean(spike_slope.subset1.slope_o,2) ; 
bw = range(spike_slope.(sdx{1}).boot_slopes(:))./50 ; 

for idx = 1:3
    hold on
    histogram(spike_slope.(sdx{1}).boot_slopes(idx,:),...
        'EdgeAlpha',0,'BinWidth',bw,'FaceColor',cm(idx,:),'FaceAlpha',0.5,...
        'Normalization','count')
    xline(meanslopes(idx),'LineWidth',2,'Color',cm(idx,:))
end

xlabel('slope')
ylabel('bootstrap count')

set(gcf,'Position',[100 100 800 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figE/' ]
mkdir(out_figdir)
filename = [out_figdir '/slope_analy_w_hist.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)
