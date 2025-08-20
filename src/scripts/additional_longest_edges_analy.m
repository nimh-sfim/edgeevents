
%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%%

% filename = [ DD.PROC '/allspk_len_indiv_' OUTSTR '.mat' ] ; 
% load(filename,'spike_len_mats')
% 
% filename = [ DD.PROC '/allspk_len_indiv80_' OUTSTR '.mat' ] ; 
% load(filename,'spike_len80_mats')
% 
% filename = [ DD.PROC '/allspk_len_indiv95_' OUTSTR '.mat' ] ; 
% load(filename,'spike_len95_mats')

subsets = {'subset1' 'subset2'} ; 
spike_len_mats = struct() ; 
spike_len80_mats= struct() ; 
spike_len95_mats= struct() ; 
spike_lenmax_mats= struct() ; 
spike_lenrange_mats = struct() ; 
spike_lenmed_mats = struct() ; 

for sdx = subsets

    spike_len_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spike_len80_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spike_len95_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spike_lenrange_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spike_lenmax_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spike_lenmed_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
       
        mm = mksq(cellfun(@mean,readdat.spike_len_cell )) ; 
        spike_len_mats.(sdx{1})(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ;

        mm2 = cellfun(@(x_) prctile(x_,[ 80 95 ]),readdat.spike_len_cell ,'UniformOutput',false) ; 

        mm80 = mksq(cellfun(@(x_) x_(1),mm2)) ; 
        mm95 = mksq(cellfun(@(x_) x_(2),mm2)) ; 

        spike_len80_mats.(sdx{1})(:,:,idx) = mm80(1:finfo.nnodes,1:finfo.nnodes) ;
        spike_len95_mats.(sdx{1})(:,:,idx) = mm95(1:finfo.nnodes,1:finfo.nnodes) ;

        ll =  mksq(cellfun(@range,readdat.spike_len_cell)) ; 
        spike_lenrange_mats.(sdx{1})(:,:,idx) = ll(1:finfo.nnodes,1:finfo.nnodes) ; 

        mm =  mksq(cellfun(@max,readdat.spike_len_cell)) ; 
        spike_lenmax_mats.(sdx{1})(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ; 

        mm =  mksq(cellfun(@median,readdat.spike_len_cell)) ; 
        spike_lenmed_mats.(sdx{1})(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ; 


    end

end


%%

subsets = {'subset1' 'subset2'} ; 
fc_conn_indiv = struct() ; 
ts_indiv = struct() ; 

for sdx = subsets

    fc_conn_indiv.(sdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        fc_conn_indiv.(sdx{1})(:,:,idx) = corr(datStr(sind).ts) ; 

    end

end

%%

addpath('/Users/faskowitzji/Documents/MATLAB/conn/')
addpath('/Users/faskowitzji/Documents/MATLAB/spm/')

%% identify the top 5% 

cmgreens = flipud(greens(100)) ; 

meanlen = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
meanlen(1:(size(meanlen,1)+1):end) = nan ;

% % dat = mean(meanlen,'omitmissing') ; 
% % parc_plot_wcolorbar(dat,surfss,annotm,...
% %     [min(dat) max(dat)],cmgreens,[100 100 600 1000])
% 
% [~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
% dat = ismember(1:200,si(1:10)) ; 
% parc_plot_wcolorbar(dat,surfss,annotm,...
%     [0 1],cmgreens,[100 100 600 1000])
% 

[~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
topdurr = ismember(1:200,si(1:10)) ; %top 5% of data

dat = double(topdurr) ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [min(dat) max(dat)],[1 1 1 ; 0.4 0.4 0.4 ],[100 100 600 1000])

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_top_nodes.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

[~,si] = sort(mean(meanlen,'omitmissing'),'ascend') ; 
bottomdurr = ismember(1:200,si(1:10)) ; %bottom 5% of data

dat = double(bottomdurr) ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [min(dat) max(dat)],[1 1 1 ; 0.6 0.6 0.6 ],[100 100 600 1000])

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_bottom_nodes.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

[~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
topdurr = ismember(1:200,si(1:10)) ; %top 5% of data

[~,si] = sort(mean(meanlen,'omitmissing'),'ascend') ; 
bottomdurr = ismember(1:200,si(1:10)) ; %bottom 5% of data

histogram(meanlen(triu(topdurr|topdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.4 .4 .4])
hold on
histogram(meanlen(triu(bottomdurr|bottomdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.6 .6 .6])
hold off
legend({'top 5% areas' 'bottom 5% areas'})

ylabel('count')
xlabel('edge event duration (edge-level)')

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_topbottom_hist.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

[~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
topdurr = ismember(1:200,si(1:40)) ; %top 20% of data

[~,si] = sort(mean(meanlen,'omitmissing'),'ascend') ; 
bottomdurr = ismember(1:200,si(1:40)) ; %bottom 20% of data

histogram(meanlen(triu(topdurr|topdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.4 .4 .4])
hold on
histogram(meanlen(triu(bottomdurr|bottomdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.6 .6 .6])
hold off
legend({'top 20% areas' 'bottom 20% areas'})

ylabel('count')
xlabel('edge event duration (edge-level)')

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_topbottom2_hist.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

[~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
topdurr = ismember(1:200,si(1:100)) ; %top 20% of data

[~,si] = sort(mean(meanlen,'omitmissing'),'ascend') ; 
bottomdurr = ismember(1:200,si(1:100)) ; %bottom 20% of data

histogram(meanlen(triu(topdurr|topdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.4 .4 .4])
hold on
histogram(meanlen(triu(bottomdurr|bottomdurr',1)),'BinWidth',0.05,'EdgeAlpha',0,...
    'FaceColor',[.6 .6 .6])
hold off
legend({'top 50% areas' 'bottom 50% areas'})

ylabel('count')
xlabel('edge event duration (edge-level)')

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_topbottom3_hist.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%%

cbig_dir = '/Users/faskowitzji/joshstuff/git_pull/CBIG/' ; 
schaefer_dir = 'stable_projects/brain_parcellation/Schaefer2018_LocalGlobal' ; 

lh_orig = [ cbig_dir schaefer_dir '/Parcellations/FreeSurfer5.3/fsaverage/surf/lh.orig' ] ;
rh_orig = [ cbig_dir schaefer_dir '/Parcellations/FreeSurfer5.3/fsaverage/surf/rh.orig' ] ;

lh_parc = [ cbig_dir schaefer_dir '/Parcellations/FreeSurfer5.3/fsaverage/label/lh.Schaefer2018_200Parcels_17Networks_order.annot' ] ; 
rh_parc = [ cbig_dir schaefer_dir '/Parcellations/FreeSurfer5.3/fsaverage/label/rh.Schaefer2018_200Parcels_17Networks_order.annot' ] ; 

[~,~,tt_L] = read_annotation(lh_parc) ; 
[~,~,tt_R] = read_annotation(rh_parc) ; 
 
cc_L = centroid_extraction_sphere(lh_orig,lh_parc) ;                        
cc_R = centroid_extraction_sphere(rh_orig,rh_parc) ;                        

rois = struct ; 

rois.sph_xyz = [ cc_L(2:end,:) ; cc_R(2:end,:) ] ; 
rois.sph_c = [ tt_L.table(2:end,1:3) ; tt_R.table(2:end,1:3) ] ./ 255 ; 

%%

cm = parula(100) ; 

[~,si] = sort(mean(meanlen,'omitmissing'),'descend') ; 
tdurr = ismember(1:200,si(1)) ; %top 5% of data

[~,si] = sort(mean(meanlen,'omitmissing'),'ascend') ; 
bdurr = ismember(1:200,si(1)) ; %top 5% of data


datrange = [min(meanlen,[],'all') max(meanlen,[],'all')] ; 

datmat{1} = tdurr|tdurr' ; 
datmat{2} = bdurr|bdurr' ; 

for ddd = 1:2

    ii = discretize(meanlen(~~datmat{ddd}),linspace(datrange(1),datrange(2),100)) ; 
    ii(isnan(ii)) = 1 ; 
    plotcols = cm(ii,:) ; 
    
    tmp = datmat{ddd} ; 
    tmp(tmp>0) = normalize(tmp(tmp>0),'range',[2 3]) ; 
    tmp(tmp==0) = 1 ; 
    rois.sph_r = datmat{ddd}' ;
    
    rois.sph_c = repmat([0.8 0.8 0.8],200,1) ; 
    
    h = viz_conn_glassbrain(datmat{ddd},plotcols ,rois) ; 
    h('bundling',7)

    out_figdir = [ './reports/figures/figR/' ] ; 
    filename = [out_figdir '/durr_edgelevel_' num2str(ddd) '.png' ] ; 
    h('print',3,filename,'-nogui')
    close(gcf)

end

% and make a colorbar
figure ; cb = colorbar ; clim(datrange)
set(gca,'Visible','off')
cb.Label.String = 'edge event duration'

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_edgelevel_cb.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

imsc_grid_comm(mean(spike_lenrange_mats.subset1,3,'omitmissing'),...
    remap_labs,1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],parc.names(g1sort))
axis square
colorbar

% cb = colorbar() ; 
% cb.Label.String = 'average event count' ; 
colormap(cmgreens)
% tmp = cb.TickLabels ; 
% cb.TickLabels = { 'non sig.' tmp{2:end} } ; 
cc = clim() ;  
clim(cc) ;
xticks([])

%% make map of 50%, 80% and 95% 

clear dat
datnames = {'50' '80' '95' } ; 
dat{1} = mean(spike_lenmed_mats.subset1,3,'omitmissing') ; 
dat{2} = mean(spike_len80_mats.subset1,3,'omitmissing') ; 
dat{3} = mean(spike_len95_mats.subset1,3,'omitmissing') ; 

for idx = 1:3

    dd = mean(dat{idx}).*finfo.TR ; 
    parc_plot_wcolorbar(dd,surfss,annotm,...
        [min(dd) max(dd)],cmgreens,[100 100 600 1000])
    
    out_figdir = [ './reports/figures/figR/' ] ; 
    mkdir(out_figdir)
    orient(gcf,'landscape')
    filename = [out_figdir '/durr_' datnames{idx} '_cortex.pdf' ] ; 
    print(filename,'-dpdf')
    close(gcf)

end

    dd = mean(mean(spike_lenrange_mats.subset1,3,'omitmissing')).*finfo.TR ; 
    parc_plot_wcolorbar(dd,surfss,annotm,...
        [min(dd) max(dd)],cmgreens,[100 100 600 1000])
    
    out_figdir = [ './reports/figures/figR/' ] ; 
    mkdir(out_figdir)
    orient(gcf,'landscape')
    filename = [out_figdir '/durr_' 'range' '_cortex.pdf' ] ; 
    print(filename,'-dpdf')
    close(gcf)


%
%%    
tiledlayout(2,1)
% 
corrmatedges = corr([tv(meanlen) tv(dat{1}) tv(dat{2}) tv(dat{3})],'type','s') ; 
corrmatmap = corr([mean(meanlen,2,'omitmissing') mean(dat{1},2,'omitmissing') ...
    mean(dat{2},2,'omitmissing') mean(dat{3},2,'omitmissing')],'type','s') ; 

% linmat = lins_ccc([tv(dat{1}) tv(dat{2}) tv(dat{3})]) ; 
% 
% % corrmat = nan(3) ; 
% % for idx = 1:3
% %     for jdx = idx+1:3
% % 
% %         corrmat(idx,jdx) = corr([tv(dat{idx}) tv(dat{jdx})]) ; 
% % 
% %     end
% % end
% 
% % tiledlayout(2,1)

nexttile()
imsc_grid_comm(corrmatedges,1:4,0.5,[1 1 1],[],{'mean' '50%' '80%' '90%'})
%clim([0.9 1])
axis square
cb = colorbar ; 
cb.Label.String = 'rank corr.' ; 
xticks('')
colormap("gray")
title('edge compare')

nexttile()
imsc_grid_comm(corrmatmap,1:4,0.5,[1 1 1],[],{'mean' '50%' '80%' '90%'})
%clim([0.9 1])
axis square
cb = colorbar ; 
cb.Label.String = 'rank corr.' ; 
xticks('')
colormap("gray")
title('map compare')



set(gcf,'Position',[100 100 800 400])

%%

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/map_corrs.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

tiledlayout(1,2)

nexttile()
imsc_grid_comm(mean(spike_lenrange_mats.subset1,3,'omitmissing').*finfo.TR,...
    remap_labs,1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],'')
axis square
colorbar

cb = colorbar() ; 
cb.Label.String = 'seconds' ; 
colormap(cmgreens) 
xticks([])
yticks([])
title('mean duration range')

nexttile()
scatter_w_rho(tv(meanlen),tv(mean(spike_lenrange_mats.subset1,3,'omitmissing')).*finfo.TR, ...
    'filled','MarkerFaceColor',[.5 .5 .5])
axis square
xlabel('mean duration')
ylabel('mean duration range')

set(gcf,'Position',[100 100 800 400])

%%

out_figdir = [ './reports/figures/figR/' ] ; 
mkdir(out_figdir)
orient(gcf,'landscape')
filename = [out_figdir '/durr_' 'range' '_matrixscatter.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


