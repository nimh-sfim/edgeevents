%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

for idx = 1:NSUBS

    disp(idx)

    filename = [DD.PROC '/' imglob '/' datStr(idx).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 

    if isfile(filename)
        disp(['already finsiehd sub: ' datStr(idx).sub ])
        continue
    end

    tmpts = datStr(idx).ts ; 
    tmpets = get_ets(tmpts) ; 

    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ; 

    mkdir([DD.PROC '/' imglob '/'])
    save(filename,...
        'spike_len_cell','spike_len_mat','-v7.3')

end

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_lengths.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
    
        spike_lengths.(sdx{1}){idx} = readdat.spike_len_cell  ; 
    
    end

end

for sdx = subsets
    tmp = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.(sdx{1}){i_}')),1:length(spike_lengths.(sdx{1}))','UniformOutput',false)') ;
    % get rid of 0 lengths
    spike_lengths.all.(sdx{1}) = nonzeros(tmp) ; 
end

%% look at spike length for each pers

% spike_prct.subset2.prct = nan(100,length(sublist.subset2)) ; 
spike_prct = struct() ;  
spike_prct.all.subset1 = prctile(spike_lengths.all.subset1, 1:100 ) ;
spike_prct.all.subset2 = prctile(spike_lengths.all.subset2, 1:100 ) ; 

for sdx = subsets

    spike_prct.(sdx{1}).prct = nan(100,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
        disp(idx)
        tmp = nonzeros(single(cell2mat(spike_lengths.(sdx{1}){idx}'))) ; 
        spike_prct.(sdx{1}).prct(:,idx) = prctile(tmp, 1:100 )  ; 
    end

end

%% 

% look at the subjects 
imagesc(spike_prct.subset1.prct(:,fliplr(sortedInd(mean(spike_prct.subset1.prct(95:100,:))))))

%% find way to naturally bin it
% look at this

plot(mean(spike_prct.subset1.prct,2),'LineWidth',2)
rr = round(mean(spike_prct.subset1.prct,2)) ; 
urr = unique(rr) ;
nums = 1:100 ; 
for idx = 1:length(unique(rr))
    hold on
    ll = rr==urr(idx) ; 
    plot(nums(ll),ones(sum(ll),1).*urr(idx),...
        Color=[0.5 0.5 0.5 0.5],LineWidth=5)
end
hold off
ylim([0 15])

%%

maxspk = 100 ; 

spike_dist = discretize(spike_lengths.all.subset1,1:1:maxspk) ; 
spike_tab = tabulate(spike_dist) ; 

% how much in lowest spike
low_bin_pct = spike_tab(1,3) ; 

% now find a high bin to get as near as possible
[~,mi] = min(abs(flipud(cumsum(flipud(spike_tab(2:end,3))))-low_bin_pct)) ; 
high_bin = mi + 1 ; % add 1 because we were aready looking at 2:end

% the low-med-high bins
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 
% 

%% 

% tiledlayout(1,2)

cm = plasma(3) ; 

clf

dd = discretize(spike_lengths.all.subset1,lowmedhigh_edges) ; 

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:40,FaceColor=cm(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    sum(dd==idx)./length(dd)
    hold on

end
hold off

xlim([1 20.5])

xlabel('spike length (sec)')
ylabel('count')

xticks((1:20)+0.5)
xticklabels(cellstr(num2str((1:20)'.*finfo.TR)))
xtickangle(45)

axis square

hist_gca = gca ;
group_lab = { ...
    strcat('short (',num2str(round(sum(dd==1)/length(dd)*100)),'%)') ... 
    strcat('intermed. (',num2str(round(sum(dd==2)/length(dd)*100)),'%)') ... 
    strcat('long (',num2str(round(sum(dd==3)/length(dd)*100)),'%)') ... 
    }

legend(hist_gca,group_lab,'Location','southeast')

aa = axes('Parent',gcf,'Position',[.48 .5 .4 .38],'Units','normalized')
box on
plot(cumsum(spike_tab(1:20,3)),'.-','MarkerSize',15,'Color',[0.5 0.5 0.5])
xlabel('spike length (sec)')
ylabel('cumulative percent')
xlim([1 20])
ylim([1 100])
axis square
xticks([1 5 10 15 20])
xticklabels(cellstr(num2str(([1 5 10 15 20])'.*finfo.TR)))

yticks([0 20 40 60 80 100])


set(gcf,'Position',[100 100 400 400])

%% 

set(gcf,'Color','w')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_hist.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% lets look at FC vs spike variability


% plot FC versus variability in lengths

% first get the variability info

[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

spike_var = struct() ; 
for sdx = subsets

    spike_var.(sdx{1}) = zeros(finfo.nnodes) ; 

    for idx = 1:length(sublist.(sdx{1}))
        disp(idx)

        dat = spike_lengths.(sdx{1}){idx}(cortmask) ; 
        tmp = mksq(cellfun(@(x_) mad(x_,0),dat)) ; 
        tmp(isnan(tmp)) = 0 ; 
            
        spike_var.(sdx{1}) = tmp + spike_var.(sdx{1}) ; 
    end

    spike_var.(sdx{1}) = spike_var.(sdx{1}) ./ length(sublist.(sdx{1})) ; 

end

%%

xdat = tv(meanfc) ; 
ydat = tv(spike_var.subset1) ; 

tiledlayout(1,2,'TileSpacing','tight')

nexttile()

imsc_grid_comm(meanfc,parc.ca(1:200),[],[],[],parc.names(1:17))
axis square
cb = colorbar ; 
cb.Label.String = 'correlation' ; 


nt2 = nexttile()

imsc_grid_comm(od_replace(spike_var.subset1,nan),parc.ca(1:200))
axis square

cb = colorbar ; 
cb.Label.String = 'spike length variability' ; 

nt2.YAxis.Visible = 'off'

set(gcf,'Position',[100 100 1000 400])

set(gcf,'Color','w')

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mean_var_matrix.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%% but also map to cortex

% MEAN FC
TL = tiledlayout(2,1)
nt1 = nexttile()

fc_degree = mean(meanfc)

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', fc_degree ,...
    'valRange',[0 max(abs(fc_degree))],...
    'cmap',parula, ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL)
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL)

hh = imagesc(fc_degree) 
cb = colorbar()
clim([0 max(abs(fc_degree))])
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cl = clim() ; 
cb.Ticks = linspace(cl(1),cl(2),5)
cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),2) , cb.TickLabels)))) ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight'

set(gcf,'Position',[100 100 600 1000])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mean_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

% FC event VAR
TL = tiledlayout(2,1)
nt1 = nexttile()

var_degree = mean(od_replace(spike_var.subset1,nan),'omitnan')

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', var_degree ,...
    'valRange',[min(abs(var_degree)) max(abs(var_degree))],...
    'cmap',parula, ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL)
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL)

hh = imagesc(var_degree) 
cb = colorbar()
clim([min(abs(var_degree)) max(abs(var_degree))])
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cl = clim() ; 
cb.Ticks = linspace(cl(1),cl(2),5)
cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),2) , cb.TickLabels)))) ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight'

set(gcf,'Position',[100 100 600 1000])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_var_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

tiledlayout(1,3,'TileSpacing','tight')


nexttile()

h = scatter_w_rho(xdat,ydat,30,'filled') ;
h.MarkerFaceColor = [0.5 0.5 0.5] ;
h.MarkerFaceAlpha = 0.25 ; 
h.MarkerEdgeAlpha = 0 ; 


axis square

xlabel('FC (correlation)')
ylabel('spike length variability')

% now colored in
nt2 = nexttile()

xbins = prctile(xdat,0:2:100) ; 

[zz,xl,mm] = norm_bin_model(xdat,xbins,ydat) ; 

h = scatter_w_rho(xdat,ydat,30,zz,'filled') ; 

cb = colorbar() ;
clim([-max(abs(zz))  max(abs(zz))])
cb.Label.String = 'within bin z-score' ;
xline(xbins) 

hold on
plot(xl,mm,'square','MarkerSize',6,'LineWidth',1,'Color','cyan')
hold off

colormap((parula))

axis square

ylabel('spike length variability')
xlabel('FC (correlation)')

% yticks([])
%nt2.YAxis.Visible = 'off'

nexttile()

h = scatter_w_rho(xdat,ydat,30,'filled') ;
h.MarkerFaceColor = [0.5 0.5 0.5] ;
h.MarkerEdgeAlpha = 0 ; 
hold on

tt = prctile(zz,95) ; 

h = scatter(xdat(zz>tt),ydat(zz>tt),30,'filled') ;
h.MarkerFaceColor = [0.2 0.2 0.2] ;
hold off

xline(xbins) 

axis square

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_zz_scatters.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% plot the relative high var based on the zz meas

TL = tiledlayout(2,1)
nt1 = nexttile()

zz_degree = mean(mksq(zz>tt)) ; 

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', zz_degree ,...
    'valRange',[min(abs(zz_degree)) max(abs(zz_degree))],...
    'cmap',[ 0.5 0.5 0.5 ; parula], ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL)
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL)

hh = imagesc(zz_degree) 
cb = colorbar()
clim([min(abs(zz_degree)) max(abs(zz_degree))])
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cl = clim() ; 
cb.Ticks = linspace(cl(1),cl(2),5)
cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),2) , cb.TickLabels)))) ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight'

set(gcf,'Position',[100 100 600 1000])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/more_spiky_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)
