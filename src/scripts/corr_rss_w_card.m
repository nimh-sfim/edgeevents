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

%% load in cardio 

loadcard = load_hcp_card('/Users/faskowitzji/joshstuff/projects/edgesort/data/processed/CARD',...
    imglob) ; 


%% now read in spike lengths to make a histogram of lengths

for sdx = subsets(1)

    corr_w_card.(sdx{1}).ccCardo = zeros(3,length(sublist.(sdx{1}))) ; 

    corr_w_card.(sdx{1}).ccardmap.short = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    corr_w_card.(sdx{1}).ccardmap.inter = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    corr_w_card.(sdx{1}).ccardmap.long = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 


    corr_w_card.(sdx{1}).ccardmapO.short = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    corr_w_card.(sdx{1}).ccardmapO.inter = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 
    corr_w_card.(sdx{1}).ccardmapO.long = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 


    corr_w_card.(sdx{1}).ccardmapTS = zeros(finfo.nnodes,length(sublist.(sdx{1}))) ; 


    corr_w_card.(sdx{1}).xcov_long_Cardo = zeros(85,length(sublist.(sdx{1}))) ; 
    corr_w_card.(sdx{1}).xcov_short_Cardo = zeros(85,length(sublist.(sdx{1}))) ; 

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


        %% card signals 
    
        cardsig = loadcard(datStr(sind).sub).values{1} ; 

        %% corr

        corr_w_card.(sdx{1}).ccCardo(:,idx) = corr(cardsig,etsrss_o) ; 

        %% corr map
        corr_w_card.(sdx{1}).ccardmapO.short(:,idx) = corr(cardsig,ets_2_ts(o1)) ; 
        corr_w_card.(sdx{1}).ccardmapO.inter(:,idx) = corr(cardsig,ets_2_ts(o2)) ; 
        corr_w_card.(sdx{1}).ccardmapO.long(:,idx) = corr(cardsig,ets_2_ts(o3)) ; 

        corr_w_card.(sdx{1}).ccardmap.short(:,idx) = corr(cardsig,ets_2_ts(dd==1)) ; 
        corr_w_card.(sdx{1}).ccardmap.inter(:,idx) = corr(cardsig,ets_2_ts(dd==2)) ; 
        corr_w_card.(sdx{1}).ccardmap.long(:,idx) = corr(cardsig,ets_2_ts(dd==3)) ; 

        corr_w_card.(sdx{1}).ccardmapTS(:,idx) = corr(cardsig,datStr(sind).ts(:,1:finfo.nnodes)) ; 

        %% xcov

        corr_w_card.(sdx{1}).xcov_short_Cardo(:,idx) = xcov(etsrss_o(:,1),cardsig,42,'normalized') ;
        corr_w_card.(sdx{1}).xcov_long_Cardo(:,idx) = xcov(etsrss_o(:,3),cardsig,42,'normalized') ;


    end
end

%%

tiledlayout(1,2)

cc = inferno(10) ; 

plotnames = {'xcov_short_Cardo' 'xcov_long_Cardo' } ; 

longerplotnames1 = {'short' 'long'} ; 
longerplotnames2 = {'cardiac' 'cardiac'} ; 

for idx = 1:2

    nexttile
    
    dat = corr_w_card.subset1.(plotnames{idx}) ; 
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

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/xcorr_w_card.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

TL = tiledlayout(2,1)
nt1 = nexttile()

dat = mean(corr_w_card.subset1.ccardmapTS,2) ; 
ccc = [ parula(100) ] ; 
vvv = [-max(abs(dat)) max(abs(dat))] ; 

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dat ,...
    'valRange',vvv,...
    'cmap',ccc, ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL)
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL)

hh = imagesc(dat) 
cb = colorbar()
clim(vvv)
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cl = clim() ; 
cb.Ticks = linspace(cl(1),cl(2),5)
%cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),5) , cb.TickLabels)))) ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight'

set(gcf,'Position',[100 100 600 1000])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figD/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_ts_corr_w_card.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

for iii = {'short' 'inter' 'long'}

    TL = tiledlayout(2,1)
    nt1 = nexttile()
    
    dat = mean(corr_w_card.subset1.ccardmapO.(iii{1}),2) ; 
    ccc = [ parula(100) ] ; 
    vvv = [-max(abs(dat)) max(abs(dat))] ; 
    
    pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dat ,...
        'valRange',vvv,...
        'cmap',ccc, ...
        'viewcMap',0,'newFig',0,'viewStr','all',...
        'parenth',TL)
    pp.Layout = nt1.Layout ; 
    
    nt2 = nexttile(TL)
    
    hh = imagesc(dat) 
    cb = colorbar()
    clim(vvv)
    hh.Visible = 'off' ;
    hh.Parent.Visible = 'off' ; 
    cb.Location = "north" ; 
    cl = clim() ; 
    cb.Ticks = linspace(cl(1),cl(2),5)
    % cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),5) , cb.TickLabels)))) ; 
    cm = colormap() ; 
    colormap(nt2,cm(2:end,:))
    
    TL.TileSpacing = 'tight'
    
    set(gcf,'Position',[100 100 600 1000])
    set(gcf,'Color','w')
    
    out_figdir = [ './reports/figures/figD/' ]
    mkdir(out_figdir)
    filename = [out_figdir '/cortex_' iii{1} '_corr_w_card.pdf' ] ; 
    print(filename,'-dpdf','-bestfit')
    close(gcf)

end

