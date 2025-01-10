%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_2.m') 

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_mats.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        spike_mats.(sdx{1}){idx} = readdat.spike_len_mat  ; 
    
    end

end

%% get counts 
 
allspike_conn_indiv = struct() ; 

for sdx = subsets

    allspike_conn_indiv.(sdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        spkmat = spike_mats.(sdx{1}){idx} ; 
        cc = count_spks(spkmat>0) ; 

        allspike_conn_indiv.(sdx{1})(:,:,idx) = mksq(cc) ; 
    end   
end

%% and get the fc for each subj

fc_conn_indiv = struct() ; 

for sdx = subsets

    fc_conn_indiv.(sdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        fc_conn_indiv.(sdx{1})(:,:,idx) = corr(datStr(sind).ts) ; 

    end

end

%% work on the figure 

ii = find(~cellfun(@isempty,regexpi(annotm('schaefer200-yeo17').combo_names,'parop'))) ; 
ind1 = ii(1) ; 

ind2 = 35 ; 
ind3 = 67 ; 

vv = zeros(200,1) ; 
vv(ind1) = 1 ; 
vv(ind2) = 2 ; 
vv(ind3) = 3 ; 

parc_plot(surfss,annotm,'schaefer200-yeo17', vv,...
        'cmap',[0.8 0.8 0.8 ; 0 0.4470 0.7410	; 0.4660 0.6740 0.1880 ; 0.8500 0.3250 0.0980  ],...
        'viewcMap',0,'newFig',0,'viewStr','lh:lat')

%%

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/example_surf.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

ts = zscore(datStr(11).ts) ; 
ee1 = ts(:,ind1).*ts(:,ind2) ; 
ee2 = ts(:,ind1).*ts(:,ind3) ; 

[c1,cc1] = count_spks(ee1,SPK_THR) ;
c2 = count_spks(ee2,SPK_THR) ;

tiledlayout(2,1)

% nexttile()
% plot(ee1,'Color',[0.8 0.8 0.8])
% ylim([-6 10]) 
% hold on
% plot(find(ee1>SPK_THR),ee1((ee1>SPK_THR)),'.')
% 
% nexttile()
% 
% plot(ee2,'Color',[0.8 0.8 0.8])
% ylim([-6 10]) 
% hold on
% plot(find(ee2>SPK_THR),ee2((ee2>SPK_THR)),'.')
% 

nexttile()
[ee1resamp,tt] = resampsig1d(ee1,1/0.72,1/0.72*10) ; 
x = tt ; 
y = ee1resamp' ; 
z = zeros(size(y)) ; 
% c = (miniets(:,idx))' ; 
c = single(ee1resamp>SPK_THR)' ; 
surface([x;x],[y;y],[z;z],[c;c],...
    'FaceColor','no','EdgeColor','interp','LineWidth',2,'EdgeAlpha',1)
colormap([0.5 0.5 0.5; 1 0 0 ])
yline(2.25,':','Color','r')

xlim([0 finfo.TR*1200])
ylim([-6 10])

ylabel('ets amp.')
xlabel('time (sec)')

text(55,8,0,['event count: ' num2str(c1)])
text(55,7,0,['correlation: ' num2str(round(mean(ee1),2))])

nexttile()
[ee2resamp,tt] = resampsig1d(ee2,1/0.72,1/0.72*10) ; 
x = tt ; 
y = ee2resamp' ; 
z = zeros(size(y)) ; 
% c = (miniets(:,idx))' ; 
c = single(ee2resamp>SPK_THR)' ; 
surface([x;x],[y;y],[z;z],[c;c],...
    'FaceColor','no','EdgeColor','interp','LineWidth',2,'EdgeAlpha',1)
colormap([0.5 0.5 0.5; 1 0 0 ])
yline(2.25,':','Color','r')

xlim([0 finfo.TR*1200])
ylim([-6 10])

ylabel('ets amp.')
xlabel('time (sec)')

text(75,8,0,['event count: ' num2str(c2)])
text(75,7,0,['correlation: ' num2str(round(mean(ee2),2))])

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/example_ets_count.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

subsets = {'subset1' 'subset2'} ; 
spike_lengths = struct() ; 

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

%%
spike_len_mats = struct() ; 

for sdx = subsets

    spike_len_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        mm = mksq(cellfun(@mean,spike_lengths.(sdx{1}){idx})) ; 
        spike_len_mats.(sdx{1})(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ;

    end

end

%% compare count vs mean

meanc = mean(allspike_conn_indiv.subset1,3) ; 
meanc = meanc(1:200,1:200) ; 

tiledlayout(1,3)

nexttile()
imsc_grid_comm(meanc,parc.ca(1:finfo.nnodes),1,[0.8 0.8 0.8],0.5)
colormap(blues)
cb = colorbar ;
cb.Label.String = "mean event count" ; 
axis square
xticks('')
yticks('')

nexttile()
imsc_grid_comm(meanfc,parc.ca(1:finfo.nnodes),1,[0.8 0.8 0.8],0.5)
colormap(blues)
cb = colorbar ;
cb.Label.String = "mean correlation" ; 
axis square
xticks('')
yticks('')

nexttile()
h = scatter_w_rho(tv(meanfc),tv(meanc),30,'filled') ;
h.MarkerFaceColor = [0.5 0.5 0.5] ;
h.MarkerFaceAlpha = 0.25 ; 
h.MarkerEdgeAlpha = 0 ; 
ylabel('count')
xlabel('correlation')


%%

set(gcf,'Position',[100 100 800 200])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/corr_vs_count.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

% %% compare count vs mean, per subject
% 
% res = zeros(length(sublist.subset1),1) ; 
% for idx = 1:length(sublist.subset1)
%     disp(idx)
%     mm = mksq(cellfun(@mean,spike_lengths.subset1{idx})) ; 
%     mm(isnan(mm)) = 0 ; 
%     ff = fc_conn_indiv.subset1(:,:,idx) ; 
%     res(idx) = corr(tv(mm),tv(ff),'type','s') ; 
% 
% end

%% compare corrvcount across cortex

res = zeros(length(sublist.subset1),finfo.nnodes) ; 
for idx = 1:length(sublist.subset1)
    disp(idx)
    mm = mksq(cellfun(@mean,spike_lengths.subset1{idx})) ; 
    mm(isnan(mm)) = 0 ; 
    ff = fc_conn_indiv.subset1(:,:,idx) ; 

    mm = mm(1:200,1:200); 
    ff = ff(1:200,1:200) ; 

    res(idx,:) = arrayfun(@(i_) corr(mm(:,i_),ff(:,i_),'type','s'),1:finfo.nnodes) ;

end

%%

parc_plot_wcolorbar(mean(res),surfss,annotm,...
    [min(mean(res)) max(mean(res))],greens(50),[100 100 600 1000])

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_corrNcount.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
%remap_colors = remaplabs(1:17,g1sort,1:17) ; 
% cmap = get_nice_yeo_cmap('grad1') ; 
cmap = greens(25) ; 
cmap = cmap(1:17,:) ; 

fcn_boxpts(mean(res)',...
    remap_labs,cmap,...
    0,parc.names(g1sort))
%([0 8.5])

set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

ylabel('mean coupling')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_corrNcount_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% finally, spike lengths

% tiledlayout(1,2)
% 
% nexttile()
% 
% ss = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
% 
% imsc_grid_comm(ss,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
% colormap(purples)
% cb = colorbar ;
% cb.Label.String = "mean event length" ; 
% axis square
% xticks('')
% yticks('')
% 
% nexttile()
% 
% %ss = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
% 
% m1 = prctile(spike_len_mats.subset1,[1],3).*finfo.TR ; 
% m2 = prctile(spike_len_mats.subset1,[99],3).*finfo.TR ; 
% 
% imsc_grid_comm(m2-m1.*finfo.TR,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
% colormap(purples)
% cb = colorbar ;
% cb.Label.String = "mean event length" ; 
% axis square
% xticks('')
% yticks('')

%%
ss = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
parc_plot_wcolorbar(mean(ss),surfss,annotm,...
    [min(mean(ss)) max(mean(ss))],viridis(50),[100 100 600 1000])

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_meanlength.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

ll = load('./data/interim/ts_autocorr.mat') ; 

autonodes =  mean(ll.acfs.subset1.map,2) ; 
autoedges = autonodes + autonodes' ; 

cbvnodes = squeeze(niftiread('./data/external/schaefer200_cbv-32k.pscalar.nii')) ; 
cbvedges = cbvnodes + cbvnodes' ; 

lll = load('./data/interim/ts_variability.mat',"vardat") ; 
varnodes = mean(lll.vardat.subset1.map,2) ; 
varedges = varnodes + varnodes' ; 

tiledlayout(2,2)

nexttile
imsc_grid_comm(ss,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
colormap(greens(100))
cb = colorbar ;
cb.Label.String = "seconds" ; 
axis square
xticks('')
yticks('')
title({ 'event length' ''} )

nexttile
imsc_grid_comm(autoedges,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
colormap(greens(100))
cb = colorbar ;
cb.Label.String = "sec.+sec." ; 
axis square
xticks('')
yticks('')
title({ 'autocorrelation' num2str(round(corr(tv(ss),tv(autoedges),'type','s'),2))} )

nexttile
imsc_grid_comm(cbvedges,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
colormap(greens(100))
cb = colorbar ;
cb.Label.String = "cbv+cbv" ; 
axis square
xticks('')
yticks('')
title({ 'cerberal blood vol.' num2str(round(corr(tv(ss),tv(cbvedges),'type','s'),2))} )

nexttile
imsc_grid_comm(varedges,parc.ca(1:finfo.nnodes),[],[0.5 0.5 0.5],0.5)
colormap(greens(100))
cb = colorbar ;
cb.Label.String = "var.+var." ; 
axis square
xticks('')
yticks('')
title({ 'time series variance' num2str(round(corr(tv(ss),tv(varedges),'type','s'),2))} )

set(gcf,'Position',[100 100 500 500])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/othermats_comp.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)