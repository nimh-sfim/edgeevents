%% start
%% TODO re-plot this

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%%

filename = [ DD.PROC '/allspk_conn_indiv_' OUTSTR '.mat' ] ; 
load(filename,'spkderiv')

%% get map of burstiness, mem

CM = reds(100) ; 

dat = mean(mean(spkderiv.subset1.burst1,3,'omitmissing')) ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [prctile(dat,0) prctile(dat,100)],flipud(CM),[100 100 600 1000])

out_figdir = [ './reports/figures/figBU/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_burst1.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

dat = mean(mean(spkderiv.subset1.burst2,3,'omitmissing')) ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [prctile(dat,0) prctile(dat,100)],flipud(CM),[100 100 600 1000])

out_figdir = [ './reports/figures/figBU/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_burst2.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


dat = mean(mean(spkderiv.subset1.mem,3,'omitmissing')) ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [prctile(dat,0) prctile(dat,100)],flipud(CM),[100 100 600 1000])

out_figdir = [ './reports/figures/figBU/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_mem.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% and system map of burst2

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
cmap = reds(25) ; 

perminds = load('./data/external/schaefer-yeo7_200node_permuted_inds.mat') ; 

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ;
dat = mean(mean(spkderiv.subset1.burst2,3,'omitmissing'))' ; 
g1 = grad_cifti(1,:)' ; 
%%

rho = corr(dat,g1,'type','s') ; 

nperm = 5e4 ; 
rng(42)
prho = nan(nperm,1) ; 
for idx = 1:nperm
    if ~mod(idx,1000) ; disp([num2str(idx/nperm*100) '%' ]) ; end
    ii = perminds.PERMS(:,randi(5e4)) ; 
    prho(idx) = corr(dat,g1(ii),'type','s') ; 
end
% spin pvalue?
pspin = (sum(abs(prho)>=abs(rho))+1)/(nperm+1) ; 

%%

tiledlayout(1,2)

nexttile()
fcn_boxpts(mean(mean(spkderiv.subset1.burst2,3,'omitmissing'))',...
    remap_labs,repmat(cmap(1,:),17,1),...
    0,parc.names(g1sort))
%([0 8.5])

ylabel('mean burstiness')
set(gca,"TickLabelInterpreter",'none')

nexttile()

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ;


[h]=scatter_w_rho(dat,grad_cifti(1,:),30,'filled')
h.MarkerFaceColor = cmap(10,:) ;     
h.MarkerEdgeColor = 'none' ;

text(0.8,0.8,['p-spin < ' num2str(round(pspin,4))],'Units','normalized')

ylabel('gradient 1 (a.u.)')
xlabel('mean burstiness')

set(gcf,'Position',[100 100 600 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figBU/' ]
mkdir(out_figdir)
filename = [out_figdir '/burst2_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

nperms = 1000 ;
NE = finfo.nnodes*(finfo.nnodes-1)/2 ; 

nullburst2 = zeros(NE,nperms) ; 

parfor idx = 1:nperms

    disp(['perm ' num2str(idx) ' of ' num2str(nperms)])

    % randomly pick cov to simulate
    randselect = randi(NSUBS); 

    ts = zscore(datStr(randselect).ts(:,1:finfo.nnodes)) ; 

    %% generate dat
    [simts_cov] = generate_phase_surrogates(ts,1,0) ; 

    tmpets = get_ets(simts_cov) ; 
    abv_thr = tmpets > SPK_THR ; 

    [~,dd] = count_spks(abv_thr) ; 
    [bb1,bb2] = spk_burst_mem(dd) ; 

    nullburst2(:,idx) = bb2 ; 

end

%% make some averages of 176

nperms = 5000 ; 
nullburstmat = zeros(finfo.nnodes,finfo.nnodes,nperms) ; 

rng(42)
for idx = 1:nperms

    disp(idx)

    ii = randsample(length(sublist.all),length(sublist.subset1),true) ; 

    bb = mean(nullburst2(:,ii),2,'omitmissing') ; 
    nullburstmat(:,:,idx) = mksq(bb) ; 

end


b2 = mean(spkderiv.subset1.burst2,3,'omitmissing') ; 

pmat = (sum(nullburstmat>=b2,3)+1)/(nperms+1) ; 

%%

sige = fdr_bh(tv(pmat),0.05) ; 
sigmask = mksq(sige) ; 
sigmask(sigmask==0) = nan ; 

sigburst = b2.*sigmask ; 

%%

dat = mean(sigburst,'omitmissing') ; 
parc_plot_wcolorbar(dat,surfss,annotm,...
    [prctile(dat,0) prctile(dat,100)],flipud(purples(50)),[100 100 600 1000])

out_figdir = [ './reports/figures/figBU/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_sigburst2.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% view top bursty edges

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

topburstthr = prctile(tv(sigburst),99) ; 
% ans =
% 
%     0.1213

wantedges = (b2>=topburstthr) & (~isnan(sigmask)) ; 

CM = reds(20) ; 

rois.sph_c = repmat([0.8 0.8 0.8],200,1) ; 

% CONN most var
h = viz_conn_glassbrain(wantedges,CM(10,:),rois) ; 

out_figdir = [ './reports/figures/figBU/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/topburst_glassbrain.png' ] ; 
h('print',3,filename,'-nogui')
close(gcf)