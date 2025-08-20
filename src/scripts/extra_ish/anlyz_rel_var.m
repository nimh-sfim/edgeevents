%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

maxspk = 100 ; 

high_bin = 4 ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

subsets = {'subset1' 'subset2'} ; 

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


%% lets save the spikevar

filename = [ DD.PROC '/spk_var_' OUTSTR '.mat' ] ; 
save(filename,'spike_var','-v7.3')

%% meas the time between spikes

spklen_names = {'short' 'inter' 'long'} ; 
spk_fano = struct() ; 

for sdx = subsets

    for ndx = 1:3
        nn = spklen_names{ndx} ; 
        spk_fano.(sdx{1}).means.(nn) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
        spk_fano.(sdx{1}).vars.(nn) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    end

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        for ndx = 1:3
            nn = spklen_names{ndx} ; 

            oo = arrayfun(@(i_) get_contact_times(~(dd(:,i_)==ndx)),1:size(dd,2),'UniformOutput',false) ; 

            spk_fano.(sdx{1}).means.(nn)(:,:,idx) = mksq(cellfun(@mean,oo)) ; 
            spk_fano.(sdx{1}).vars.(nn)(:,:,idx) = mksq(cellfun(@var,oo)) ; 
        end

    end
end

%%

filename = [ DD.PROC '/spk_fano_' OUTSTR '.mat' ] ; 
save(filename,'spk_fano','-v7.3')

%% plot fano stuff?



%% FC vs spike length variability

tiledlayout(1,2,'TileSpacing','tight')

nt1 = nexttile()

imsc_grid_comm(meanfc,parc.ca(1:200),[],[],1,parc.names(1:17))
axis square
cb = colorbar ; 
cb.Label.String = 'correlation' ; 
xticks([])
colormap(nt1,gray(100))

nt2 = nexttile() ; 

imsc_grid_comm(od_replace(spike_var.subset1,nan),parc.ca(1:200),[],[],1)
axis square

cb = colorbar ; 
cb.Label.String = 'spike length variability' ; 

nt2.YAxis.Visible = 'off'
xticks([])

colormap(nt2,gray(100))

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mean_var_matrix.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%% map fc and event var to cortex

fc_degree = mean(meanfc) ;

parc_plot_wcolorbar(fc_degree,surfss,annotm,...
    [0 max(abs(fc_degree))],fblue(100),[100 100 600 1000])

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mean_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

var_degree = mean(od_replace(spike_var.subset1,nan),'omitnan') ; 

parc_plot_wcolorbar(var_degree,surfss,annotm,...
    [0 max(abs(var_degree))],fblue(100),[100 100 600 1000])

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_var_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% MAKE ZZ -- residualize variability per FC

xdat = tv(meanfc) ; 
ydat = tv(spike_var.subset1) ; 

xbins = prctile(xdat,0:5:100) ; % just use simple percentiles, all bins have
                                % same numeber of edges

[eVar,xbinmid,binmean] = norm_bin_model(xdat,xbins,ydat) ; 

%% lets to a lil testing on the plot

%grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 
perminds = load('./data/external/schaefer-yeo7_200node_permuted_inds.mat') ; 

eVarMat = mksq(eVar) ; 
bb = get_blocky(eVarMat,parc.ca(1:200)) ; 
nperm = 1e4 ; 
nblocks = size(bb,1) ; 
permbb = nan(nblocks,nblocks,nperm) ; 
rng(42)
for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perminds.PERMS(:,randi(5e4)) ; 
    permbb(:,:,idx) = get_blocky(eVarMat(ii,ii),parc.ca(1:200)) ; 
end

phigh = sum(bb<=permbb,3)./nperm ; 
plow = sum(bb>=permbb,3)./nperm ;

%%

mmask = logical(triu(ones(17),0)) ; 

tmp = fdr_bh(phigh(mmask),0.05) ; 
sighigh_eVar_blocks = zeros(17) ; 
sighigh_eVar_blocks(mmask) = tmp ; 

tmp = fdr_bh(plow(mmask),0.05) ; 
siglow_eVar_blocks = zeros(17) ; 
siglow_eVar_blocks(mmask) = tmp ; 


%%

tiledlayout(1,3,'TileSpacing','tight')

CM = thermal(100) ; 

nexttile()

h = scatter_w_rho(xdat,ydat,30,'filled') ;
h.MarkerFaceColor = [0.5 0.5 0.5] ;
h.MarkerFaceAlpha = 0.25 ; 
h.MarkerEdgeAlpha = 0 ; 

axis square

xlabel('FC (correlation)')
ylabel('spike length variability')

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

ylabel('spike length variability')
xlabel('FC (correlation)')

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

xline(xbins) 
ylabel('spike length variability')

axis square
set(gcf,'Position',[100 100 1200 600])
set(gcf,'Color','w')
orient(gcf,'landscape')

%%

out_figdir = [ './reports/figures/figB/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_zz_scatters.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% plot matrix of zscore var vals + yeo communities 

TL = tiledlayout(1,3,'TileSpacing','tight')

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

nexttile(TL)

imsc_grid_comm(mksq(eVar),remap_labs,[],[1 1 1],1,parc.names(1:17))
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

out_figdir = [ './reports/figures/figB/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/spike_zz_mat.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% sys

tiledlayout(1,2)

nexttile()

g1sort = [ 13 17 14 8 16 11 7 15 12 10 1 3 6 9 2 4 5 ] ; 

remap_nodes = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
[datplot_ca,datplot_sinds] = sort(remap_nodes) ;  

cmapremap = remaplabs(1:17,g1sort,1:17) ; 
%cmap = get_nice_yeo_cmap('grad1') ; 
%cmap = cmap(cmapremap,:) ; 
cmap = viridis(17) ; 

% get the within-community edge vals
[uu,vv] = find(triu(ones(finfo.nnodes),1)) ; 
p_uu = remap_nodes(uu) ; p_vv = remap_nodes(vv) ;

vals = cell(17,1) ; 
vallabs = cell(17,1) ; 
for idx = 1:17
    vals{idx} = eVar((p_uu==idx) & (p_vv==idx)) ; 
    vallabs{idx} = ones(size(vals{idx})).*idx ; 
end

rng(42)
fcn_boxpts(cell2mat(vals),cell2mat(vallabs),cmap,0,parc.names(g1sort))

xtickangle(45)
set(gca,'TickLabelInterpreter','none')

axis square

title('within-system edge values')
ylabel('within bin z-score')


nexttile()

[aa,bb,cc] = perm_mat_anova1(mksq(eVar),remap_nodes,1000) ; 

mm = multcompare(cc,'Alpha',0.001,'Display','off') ; 
%mmm = mksq(mm(:,6)) ; 
trilm = logical(tril(ones(17),-1)) ;
mmm = zeros(17) ; 
mmm(trilm) = mm(:,6) ;
mmm = mmm + mmm' ; 

h = imsc_grid_comm(mmm,1:17,0.5,[1 1 1],[1 1 1],parc.names(g1sort)) ; 
axis square
sigthr = 1/(17*16/2) ; 
tmp = (mmm <= sigthr)+.5 ; 
tmp(~~eye(size(tmp))) = 0 ; 
h.AlphaData = tmp ;

cb = colorbar() ; 
cb.Label.String = 'p-value' ; 

title("Tukey's HSD")


set(gcf,'Position',[100 100 600 600])
set(gcf,'Color','w')
orient(gcf,'landscape')


%%

out_figdir = [ './reports/figures/figB/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/spike_zz_sys.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% plot the relative high var based on the zz meas

mycolmap  = thermal(200)  ; 

zz_degree = mean(mksq(eVar>=topvarthr)) ; 

parc_plot_wcolorbar(zz_degree,surfss,annotm,...
    [min(abs(zz_degree)) max(abs(zz_degree))], ...
    [ 0.5 0.5 0.5 ; mycolmap(size(mycolmap,1)/2:end,:)],[100 100 600 1000])

out_figdir = [ './reports/figures/figB/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/yesvar_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

zz_degree = mean(mksq(eVar<=bottomvarthr)) ; 

parc_plot_wcolorbar(zz_degree,surfss,annotm,...
    [min(abs(zz_degree)) max(abs(zz_degree))], ...
    [ 0.5 0.5 0.5 ; flipud(mycolmap(1:size(mycolmap,1)/2,:))],[100 100 600 1000])

out_figdir = [ './reports/figures/figB/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/novar_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
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

% lh_parc2 = '/Users/faskowitzji/joshstuff/git_pull/parc_plotter/data/fsaverage/label/lh.schaefer200-yeo17.annot' ; 
% rh_parc2 = '/Users/faskowitzji/joshstuff/git_pull/parc_plotter/data/fsaverage/label/rh.schaefer200-yeo17.annot' ; 
% 
% [~,~,tt_L2] = read_annotation(lh_parc2) ; 
% [~,~,tt_R2] = read_annotation(rh_parc2) ; 
% 
% cc_L2 = centroid_extraction_sphere(lh_orig,lh_parc2) ;                        
% cc_R2 = centroid_extraction_sphere(rh_orig,rh_parc2) ;                        
% 
% rois2 = struct ; 
% 
% rois2.sph_xyz = [ cc_L2(2:end,:) ; cc_R2(2:end,:) ] ; 
% rois2.sph_c = [ tt_L2.table(2:end,1:3) ; tt_R2.table(2:end,1:3) ] ./ 255 ; 
% 

%%

CM = thermal(100) ; 

eVarMat = mksq(eVar) ;
mostvar_ind = eVar>=prctile(eVar,99)  ; 
mostvar_sqr = triu(mksq(mostvar_ind),1) ; 

leastvar_ind = eVar<=prctile(eVar,1)  ; 
leastvar_sqr = triu(mksq(leastvar_ind),1) ; 

% cm = parula(100) ; 
% cm = cm(51:100,:) ; 
% assign all the values to colormap ind
% mostvar_col =  cm(discretize(zz(mostvar_ind),50),:) ; 

rois.sph_c = repmat([0.8 0.8 0.8],200,1) ; 

tmp = sum(mostvar_sqr+mostvar_sqr') ; 
tmp(tmp>0) = normalize(tmp(tmp>0),'range',[2 5]) ; 
tmp(tmp==0) = 1 ; 
rois.sph_r = tmp' ; 

% CONN most var
h = viz_conn_glassbrain(mostvar_sqr+mostvar_sqr',CM( floor(size(CM,1)*0.8),:) ,rois) ; 

out_figdir = [ './reports/figures/figB/' ] ; 
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

out_figdir = [ './reports/figures/figB/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/novar_glassbrain.png' ] ; 
h('print',3,filename,'-nogui')
close(gcf)

% %pp = parula(10) ; 
% 
% eVarMat = mksq(eVar) ;
% mostvar_ind = eVar>=prctile(eVar,99)  ; 
% mostvar_sqr = triu(mksq(mostvar_ind),1) ; 
% 
% cm = parula(100) ; 
% %cm = cm(51:100,:) ; 
% 
% fcvec = tv(meanfc) ; 
% fccol = cm(discretize(fcvec,linspace(min(fcvec),max(fcvec),100)),:) ; 
% 
% scatter(fcvec(mostvar_ind),fcvec(mostvar_ind),40,fccol(mostvar_ind,:),'filled')
% 
% 
% scatter(1:199,fcvec(mostvar_ind),40,fccol(mostvar_ind,:),'filled')
% 
% % assign all the values to colormap ind
% % mostvar_col =  cm(discretize(fcvec(mostvar_ind),50),:) ; 
% 
% rois.sph_c = repmat([0.8 0.8 0.8],200,1) ; 
% 
% sum(mostvar_ind)
% 
% h = viz_conn_glassbrain(mostvar_sqr,fccol(mostvar_ind,:) ,rois) ; 

%% find significant blocks

zzsig = struct() ; 
zznames = {'mostvar' 'leastvar' };
zzdat = cell(2,1) ; 
zzdat{1} = (mostvar_sqr+mostvar_sqr') ; 
zzdat{2} = (leastvar_sqr+leastvar_sqr') ;

for sdx = 1:2

    mat = zzdat{sdx} ; 
    bb = get_blocky(mat,parc.ca(1:200)) ; 
    
    nperm = 1e4 ; 
    nblocks = size(bb,1) ; 
    permbb = nan(nblocks,nblocks,nperm) ; 
    rng(42)
    for idx = 1:nperm
        if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
        ii = perminds.PERMS(:,randi(5e4)) ; 
        permbb(:,:,idx) = get_blocky(mat(ii,ii),parc.ca(1:200)) ; 
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
%     {'ContA'    }    {'VisCent'}
%     {'DorsAttnA'}    {'VisCent'}

%% also boxplot for the mostvar edges

ss = sum(mostvar_sqr + mostvar_sqr') ; 
ca = parc.ca(1:200) ; 

%% 

binedges = [ [xbins(17) xbins(18)]; ...
    [xbins(18) xbins(19)] ; ...
    [xbins(19) xbins(20)] ] ;
    
fcvec = tv(meanfc) ; 

for idx = 1:3

    be = binedges(idx,:) ; 

    edgesinbin = fcvec>=be(1) & fcvec<be(2) ; 
    
    eVarColors = CM(discretize(eVar,linspace(min(eVar),max(eVar),100)),:) ; 
    
    tmp = sum(mksq(edgesinbin)) ; 
    tmp(tmp>0) = normalize(tmp(tmp>0),'range',[2 3]) ; 
    tmp(tmp==0) = 1 ; 
    rois.sph_r = tmp' ;
    
    h = viz_conn_glassbrain(triu(mksq(edgesinbin),1),eVarColors(edgesinbin,:) ,rois) ; 
    
    out_figdir = [ './reports/figures/figB/' ] ; 
    mkdir(out_figdir)
    filename = [out_figdir '/bin_' num2str(idx) '_glassbrain.png' ] ; 
    h('print',3,filename,'-nogui')
    close(gcf)

end

%% nuthin here

% binedges = [ 
%     [xbins(15) xbins(16)] ; ...
%     [xbins(16) xbins(17)] ; ...
%     [xbins(17) xbins(18)] ; ...
%     [xbins(18) xbins(19)] ; ...
%     [xbins(19) xbins(20)] ] ;
% 
% fcvec = tv(meanfc) ; 
% 
% % hh = zeros(finfo.nnodes,1) ; 
% % hh(101:200) = 1 ; 
% hh = zeros(finfo.nnodes,finfo.nnodes) ; 
% hh(1:100,1:100) = 1 ; 
% hh(101:200,101:200) = 1 ; 
% 
% for idx = 1:5
% 
%     be = binedges(idx,:) ; 
% 
%     edgesinbin = fcvec>=be(1) & fcvec<be(2) ; 
% 
%     eVarColors = CM(discretize(eVar,linspace(min(eVar),max(eVar),100)),:) ; 
% 
%     tmp = (mksq(edgesinbin.*eVar)) ; 
% 
%     % imsc_grid_comm(tmp,parc.ca(1:200))
%     % title(idx)
%     % axis square
% 
%     % parc_plot(surfss,annotm,'schaefer200-yeo17',...
%     %     sum(tmp),'valRange',[-30 30],...
%     %     'cmap',parula(100),'viewcMap',0,'newFig',0,'viewStr','all')    
%     % waitforbuttonpress
% 
%     %fcn_boxpts([sum(tmp.*(hh==1)) sum(tmp.*(hh==0))]',[ ones(200,1) ; ones(200,1).*2])
%     parallelcoords([sum(tmp.*(hh==1)) ; sum(tmp.*(hh==0))]')
%     title(idx)
%     waitforbuttonpress
%     clf
% 
% end
% 
