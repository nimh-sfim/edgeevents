%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%% load up the data you need

filename = [ DD.PROC '/allspk_conn_indiv_' OUTSTR '.mat' ] ; 
load(filename)

filename = [ DD.PROC '/allspk_len_indiv_' OUTSTR '.mat' ] ; 
load(filename)

%%

meanc = mean(allspike_conn_indiv.subset1,3) ; 
meanc = meanc(1:200,1:200) ; 

meanlen = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
meanlen(1:size(meanlen,1)+1:end) = nan ;

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

%%

CM = flipud(rdbu(100)) ; 

xdat = tv(meanc) ; 
ydat = tv(meanlen) ; 

xbins = prctile(xdat,0:4:100) ; % just use simple percentiles, all bins have
                               % same numeber of edges
[eVar,xbinmid,binmean] = norm_bin_model(xdat,xbins,ydat) ; 

% h = scatter(tv(meanc),tv(meanlen),30,eVar,'filled') ; 
% 
% xlabel('mean count')
% ylabel('mean dur.')
% 
% cb = colorbar() ;
% minmaxlim = prctile(abs(eVar),99) ; 
% clim([-minmaxlim  minmaxlim])
% cb.Label.String = 'within bin z-score' ;
% xline(xbins) 
% 
% hold on
% plot(xbinmid,binmean,'square','MarkerSize',6,'LineWidth',1,'Color','cyan')
% hold off
% 
% axis square
% 
% colormap(CM)

%%
%%

perminds = load('./data/external/schaefer-yeo7_200node_permuted_inds.mat') ; 

eVarMat = mksq(eVar) ; 
bb = get_blocky(eVarMat,remap_labs) ; 
nperm = 1e4 ; 
nblocks = size(bb,1) ; 
permbb = nan(nblocks,nblocks,nperm) ; 
rng(42)
for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perminds.PERMS(:,randi(5e4)) ; 
    permbb(:,:,idx) = get_blocky(eVarMat(ii,ii),remap_labs) ; 
end

phigh = sum(bb<=permbb,3)./nperm ; 
plow = sum(bb>=permbb,3)./nperm ;

mmask = logical(triu(ones(17),0)) ; 

tmp = fdr_bh(phigh(mmask),0.05) ; 
sighigh_eVar_blocks = zeros(17) ; 
sighigh_eVar_blocks(mmask) = tmp ; 

tmp = fdr_bh(plow(mmask),0.05) ; 
siglow_eVar_blocks = zeros(17) ; 
siglow_eVar_blocks(mmask) = tmp ; 

%%

tiledlayout(1,3,'TileSpacing','tight')

CM = rdbu(100) ; 

nexttile()

h = scatter_w_rho(xdat,ydat,30,'filled') ;
h.MarkerFaceColor = [0.5 0.5 0.5] ;
h.MarkerFaceAlpha = 0.25 ; 
h.MarkerEdgeAlpha = 0 ; 

axis square

ylabel('mean even duration')
xlabel('mean event count')

% now colored in
nexttile()

h = scatter(xdat,ydat,30,eVar,'filled') ; 

cb = colorbar() ;
clim([-max(abs(eVar))  max(abs(eVar))])
cb.Label.String = 'within bin z-score' ;
xline(xbins) 

hold on
plot(xbinmid,binmean,'square','MarkerSize',6,'LineWidth',1,'Color','cyan')
hold off

colormap(CM)

axis square

ylabel('mean even duration')
xlabel('mean event count')

% yticks([])
%nt2.YAxis.Visible = 'off'

nexttile()

h = scatter(xdat,ydat,30,'filled') ;
h.MarkerFaceColor = [0.8 0.8 0.8] ;
h.MarkerEdgeAlpha = 0 ; 
hold on

topvarthr = prctile(eVar,99) ; 
bottomvarthr = prctile(eVar,1) ; 

h = scatter(xdat(eVar>=topvarthr),ydat(eVar>=topvarthr),30,'filled') ;
% h.MarkerFaceColor = [0.2 0.2 0.2] ;
h.MarkerFaceColor = CM( floor(size(CM,1)*0.8),:) ;

h = scatter(xdat(eVar<=bottomvarthr),ydat(eVar<=bottomvarthr),30,'filled') ;
% h.MarkerFaceColor = [0.2 0.2 0.2] ;
h.MarkerFaceColor = CM( floor(size(CM,1)*0.2),:) ;

hold off

%xline(xbins
 
ylabel('mean even duration')
xlabel('mean event count')

axis square
set(gcf,'Position',[100 100 1200 600])
set(gcf,'Color','w')
orient(gcf,'landscape')

%%

out_figdir = [ './reports/figures/figCnC/' ]
mkdir(out_figdir)
filename = [out_figdir '/cnc_eVar_scatters.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% plot matrix of zscore var vals + yeo communities 

TL = tiledlayout(1,3,'TileSpacing','tight')

nexttile(TL)

imsc_grid_comm(mksq(eVar),remap_labs,[],[1 1 1],1,parc.names(g1sort))
axis square
colormap(CM)

cb = colorbar() ;
clim([-max(abs(eVar))  max(abs(eVar))])
cb.Label.String = 'within bin z-score' ;

xticks([])

nexttile(TL)

imsc_grid_comm(get_blocky(mksq(eVar),remap_labs),1:17,[],[1 1 1],1,[])
axis square
colormap(CM)
clim([-max(abs(eVar))  max(abs(eVar))])
cb = colorbar() ;
clim([-max(abs(eVar))  max(abs(eVar))])
cb.Label.String = 'within bin z-score' ;

xticks([])

nexttile(TL)

h = imsc_grid_comm(get_blocky(mksq(eVar) ,remap_labs),...
    1:17,0.5,[.8 .8 .8],[0.8 .8 .8],[]) ;
h.AlphaData = sighigh_eVar_blocks + siglow_eVar_blocks' ; 
axis square
colormap(CM)
clim([-max(abs(eVar))  max(abs(eVar))])
% cb = colorbar() ;
% clim([-max(abs(eVar))  max(abs(eVar))])
% cb.Label.String = 'within bin z-score' ;

[uu,vv] = find(sighigh_eVar_blocks) ; 
for idx = 1:length(uu)
    hold on
    imsc_addsquare(1:17==uu(idx),1:17==vv(idx),0,CM( floor(size(CM,1)*0.8),:))
end
hold off

[uu,vv] = find(siglow_eVar_blocks') ; 
for idx = 1:length(uu)
    hold on
    imsc_addsquare(1:17==uu(idx),1:17==vv(idx),0,CM( floor(size(CM,1)*0.2),:))
end
hold off

xticks([])
yticks([])

set(gcf,'Position',[100 100 1400 1000])
set(gcf,'Color','w')
orient(gcf,'landscape')

%% 

out_figdir = [ './reports/figures/figCnC/' ]
mkdir(out_figdir)
filename = [out_figdir '/cnc_eVar_mats.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

addpath('/Users/faskowitzji/Documents/MATLAB/spm12/')
addpath('/Users/faskowitzji/Documents/MATLAB/conn/')

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

topvarthr = prctile(eVar,98) ; 
bottomvarthr = prctile(eVar,2) ; 

eVarMat = mksq(eVar) ;
mostvar_ind = eVar>=topvarthr  ; 
mostvar_sqr = triu(mksq(mostvar_ind),1) ; 

leastvar_ind = eVar<=bottomvarthr  ; 
leastvar_sqr = triu(mksq(leastvar_ind),1) ; 

rois.sph_c = repmat([0.8 0.8 0.8],200,1) ; 

tmp = sum(mostvar_sqr+mostvar_sqr') ; 
tmp(tmp>0) = normalize(tmp(tmp>0),'range',[2 5]) ; 
tmp(tmp==0) = 1 ; 
rois.sph_r = tmp' ; 

mm = logical(mostvar_sqr+mostvar_sqr') ; 
mm(~~eye(size(mm))) = 0 ;

% CONN most var
h = viz_conn_glassbrain(mm,CM( floor(size(CM,1)*0.8),:) ,rois) ; 

out_figdir = [ './reports/figures/figCNC/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/yesvar_glassbrain.png' ] ; 
h('print',3,filename,'-nogui')
close(gcf)


tmp = sum(leastvar_sqr+leastvar_sqr') ; 
tmp(tmp>0) = normalize(tmp(tmp>0),'range',[2 5]) ; 
tmp(tmp==0) = 1 ; 
rois.sph_r = tmp' ; 

% CONN least var
h = viz_conn_glassbrain(leastvar_sqr+leastvar_sqr',CM( floor(size(CM,1)*0.2),:) ,rois) ; 

out_figdir = [ './reports/figures/figCnC/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/novar_glassbrain.png' ] ; 
h('print',3,filename,'-nogui')
close(gcf)

%%

zz_degree = mean(mostvar_sqr+mostvar_sqr') ; 
zz_colmap = [ [0.5 0.5 0.5] ; interp_cmap_lin([0.8 0.8 0.8],CM( floor(size(CM,1)*0.8),:),100) ] ; 

parc_plot_wcolorbar(zz_degree,surfss,annotm,...
    [min(abs(zz_degree)) max(abs(zz_degree))], ...
    zz_colmap,[100 100 600 1000])

out_figdir = [ './reports/figures/figCnC/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/yesvar_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

zz_degree = mean(leastvar_sqr+leastvar_sqr') ; 
zz_colmap = [ [0.5 0.5 0.5] ; interp_cmap_lin([0.8 0.8 0.8],CM( floor(size(CM,1)*0.2),:),100) ] ; 

parc_plot_wcolorbar(zz_degree,surfss,annotm,...
    [min(abs(zz_degree)) max(abs(zz_degree))], ...
    zz_colmap,[100 100 600 1000])

out_figdir = [ './reports/figures/figCnC/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/novar_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% find significant blocks

zzsig = struct() ; 
zznames = {'mostvar' 'leastvar' };
zzdat = cell(2,1) ; 
zzdat{1} = (mostvar_sqr+mostvar_sqr') ; 
zzdat{2} = (leastvar_sqr+leastvar_sqr') ;

for sdx = 1:2

    mat = zzdat{sdx} ; 
    bb = get_blocky(mat,remap_labs) ; 
    
    nperm = 1e4 ; 
    nblocks = size(bb,1) ; 
    permbb = nan(nblocks,nblocks,nperm) ; 
    rng(42)
    for idx = 1:nperm
        if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
        ii = perminds.PERMS(:,randi(5e4)) ; 
        permbb(:,:,idx) = get_blocky(mat(ii,ii),remap_labs) ; 
    end
    
    phigh = sum(bb<=permbb,3)./nperm ; 
    plow = sum(bb>=permbb,3)./nperm ;
    
    mmask = logical(triu(ones(17),0)) ; 
    
    tmp = fdr_bh(phigh(mmask),0.05) ; 
    tmpblocky = zeros(17) ; 
    tmpblocky(mmask) = tmp ; 
    
    zzsig.(zznames{sdx}).high = tmpblocky ;

    tmp = fdr_bh(plow(mmask),0.05) ; 
    tmpblocky = zeros(17) ; 
    tmpblocky(mmask) = tmp ; 
    
    zzsig.(zznames{sdx}).low = tmpblocky ;

end

[uu,vv] = find(zzsig.mostvar.high) ; 
disp('sig connections in mostvar:')
disp([ parc.names(uu) parc.names(vv) ] )
% sig connections in mostvar:
%     {'DefaultB'  }    {'DorsAttnA'   }
%     {'Limbic_OFC'}    {'SalVentAttnA'}

%% what would the 'regular' top and bottom length egdes look like

toplenedges = tv(meanlen)>=prctile(tv(meanlen),98) ; 
bottomlenedges = tv(meanlen)<=prctile(tv(meanlen),2) ; 


% CONN most var
h = viz_conn_glassbrain(mksq(toplenedges),CM( floor(size(CM,1)*0.8),:) ,rois) ; 


zz_degree = mean(mksq(toplenedges)) ; 
zz_colmap = [ [0.5 0.5 0.5] ; interp_cmap_lin([0.8 0.8 0.8],CM( floor(size(CM,1)*0.8),:),100) ] ; 

parc_plot_wcolorbar(zz_degree,surfss,annotm,...
    [min(abs(zz_degree)) max(abs(zz_degree))], ...
    zz_colmap,[100 100 600 1000])

zz_degree2 = mean(mksq(eVar>topvarthr)) ; 
zz_colmap = [ [0.5 0.5 0.5] ; interp_cmap_lin([0.8 0.8 0.8],CM( floor(size(CM,1)*0.8),:),100) ] ; 

parc_plot_wcolorbar(zz_degree2,surfss,annotm,...
    [min(abs(zz_degree2)) max(abs(zz_degree2))], ...
    zz_colmap,[100 100 600 1000])

%%

%% what if you look at event length, outside of on-diagonal 

newmat = zeros(finfo.nnodes,finfo.nnodes) ; 
ca = parc.ca(1:finfo.nnodes) ; 


newmat = meanlen .* 1 ; 
newmat(ca==ca') = nan ; 

dat = mean(newmat,'omitnan') ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [min(dat) max(dat)],flipud(greens(100)),[100 100 600 1000])

dat = mean(meanlen,'omitnan') ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [min(dat) max(dat)],flipud(greens(100)),[100 100 600 1000])

