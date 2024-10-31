%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

%% load the spike conn

filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
load(filename)

%% look at long vs short scatter again, and extract just the super long edges

% dat1 = tv(spike_conn.subset1.short) ; 
% dat2 = tv(spike_conn.subset1.long) ; 
% 
% tld = arrayfun(@(i_) point_to_line([dat1(i_) dat2(i_) 0],[0 0 0],[1 1 0]),1:length(dat1)) ; 
% 
% longest_edges = tld>3 ; 


%% do some stat testing

% load up the nulls
filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
load(filename)

%% boots

nboot = 1000 ; 

nsubs = length(sublist.all) ; 
subsetsz = length(sublist.subset1) ; 

nperms = length(simmat.keepcov) ; 
rng(42)

surrA = struct() ;
surrB = struct() ;
lennames = {'short' 'inter' 'long'} ; 

% and record for this boot
for sdx = 1:3
    surrA.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot); 
    surrB.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot) ; 
end

for idx = 1:nboot
    
    disp(idx)

    % pickout the surr data to average
    bootinds = randsample(1:nperms,subsetsz,true) ; 

    tmpA = zeros(3,finfo.nnodes*(finfo.nnodes-1)/2) ; 
    tmpB = zeros(3,finfo.nnodes*(finfo.nnodes-1)/2)  ; 

    for bdx = bootinds 

        aa = arrayfun(@(i_)simmat.nocov{bdx}(i_,:),1:3,'UniformOutput',false) ; 
        bb = arrayfun(@(i_)simmat.keepcov{bdx}(i_,:),1:3,'UniformOutput',false) ; 

        for sdx = 1:3
            tmpA(sdx,:) = tmpA(sdx,:) + aa{sdx}  ; 
            tmpB(sdx,:) = tmpB(sdx,:) + bb{sdx}  ; 
        end

    end
    
    % make mean
    for sdx = 1:3
        tmpA(sdx,:) = tmpA(sdx,:) ./ subsetsz  ; 
        tmpB(sdx,:) = tmpB(sdx,:) ./ subsetsz  ; 
    end

    % and record for this boot
    for sdx = 1:3
        surrA.(lennames{sdx})(:,:,idx) = mksq(tmpA(sdx,:)) ; 
        surrB.(lennames{sdx})(:,:,idx) = mksq(tmpB(sdx,:)) ; 
    end

end

%%

sigmoreA = struct() ; 
sigmoreB = struct() ; 

siglessA = struct() ; 
siglessB = struct() ; 

for sdx = 1:3

    sigmoreA.(lennames{sdx}) = spike_conn.subset1.(lennames{sdx}) >= ...
        prctile(surrA.(lennames{sdx}),99.9,3) ; 
    sigmoreB.(lennames{sdx}) = spike_conn.subset1.(lennames{sdx}) >= ...
        prctile(surrB.(lennames{sdx}),99.9,3) ; 

    siglessA.(lennames{sdx}) = spike_conn.subset1.(lennames{sdx}) <= ...
        prctile(surrA.(lennames{sdx}),0.1,3) ; 
    siglessB.(lennames{sdx}) = spike_conn.subset1.(lennames{sdx}) <= ...
        prctile(surrB.(lennames{sdx}),0.1,3) ; 

end


%% sig test 

CM = flipud(tempo(100)) ; 

TL = tiledlayout(1,3)

nt = nexttile()

% from a different approach
%
% longedge_thr = min(spike_conn.subset1.long(sig_longedges)) ; % prctile(dat,95) ; 
% ll = linspace(0,longedge_thr,20) ; 
% ss = ll(2) ; 
% 
% hh = histogram(dat,'BinEdges', ll ,FaceColor=[0.8 0.8 0.8],...
%     Normalization='count')
% hold on
% histogram(dat(dat>=longedge_thr),'FaceColor',[0.2 0.2 0.2],...
%     'BinEdges',ll(end):ss:max(dat),...
%     Normalization='count')
% 
% ylabel('count')
% xlabel('avg. num. long')
% xx = xticks() ; 
% xticklabels(cellstr(num2str(round(xx.*finfo.TR,2)')))
% 
% line([ longedge_thr longedge_thr ], [0 1750] ,'Color','red')
% text(longedge_thr,2000, ...
%     [ { 'surrogate-based' ['thr: ' num2str(round(longedge_thr*finfo.TR,2)) ] } ],...
%     'HorizontalAlignment','center')
% 
% axis square

% aa = axes('Parent',gcf,'Position',[.48 .5 .4 .38],'Units','normalized')
% ss = scatter((dat1),(dat2),10,tld>3,"filled") ; 
% cc = inferno(10) ; 
% colormap(nt,[0.8 0.8 0.8 ; cc(2,:) ])
% rl = refline(1,0) ; 
% rl.LineWidth = 2 ;
% colormap(inferno(2))

imsc_grid_comm(sigmoreA.long.*spike_conn.subset1.long,parc.ca(1:200),1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],parc.names(1:17))
axis square

xticks([])
cb = colorbar() ; 
cb.Label.String = 'average event count' ; 
colormap([ 1 1 1 ; CM])
tmp = cb.TickLabels ; 
cb.TickLabels = { 'non sig.' tmp{2:end} } ; 
cc = clim() ;  
clim(cc) ;

nt2 = nexttile() ;

%longest_edges = dat >= longedge_thr ; 

bb = get_blocky(sigmoreA.long,parc.ca(1:200))
imsc_grid_comm(bb,1:17,...
    1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],[])
set(gca,'TickLength',[ 0 0])
cb = colorbar()
cb.Label.String = 'density of sig. edges' ; 
xticks([]) 
yticks([])

% xlabel('btwn system density')

axis square

colormap(nt2,[ 1 1 1 ; flipud(tempo(100))])

nexttile()

histogram(nonzeros(tv(sigmoreA.long.*spike_conn.subset1.long)),...
    'Normalization','probability','FaceColor',CM(80,:),'EdgeAlpha',0)
hold on
histogram(nonzeros(tv(~sigmoreA.long.*spike_conn.subset1.long)),...
    'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0)

legend({'sig. edges' 'non-sig. edges'},'Location','northeastoutside')

xlabel('average event length')
ylabel('prob.')

axis square

set(gcf,'Position',[100 100 1000 400])

%% 

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_spike_hist.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%% glassbrain
 
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

% too many edges for plot
% rois.sph_r = normalize(sum(sigmoreA.long), 'range',[1,3]) ; 
% viz_conn_glassbrain(sigmoreA.long,[ 0.5 .5 .5],rois)
% 
% out_figdir = [ './reports/figures/figC/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/longest_spike_glassbrain.png' ] ; 
% h('print',3,filename,'-nogui')
% close(gcf)

parc_plot_wcolorbar(sum(sigmoreA.long),surfss,annotm,...
    [0 max(abs(sum(sigmoreA.long)))],CM,[100 100 600 1000])

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longsig_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)


% 
% set(gcf,'Position',[100 100 600 400])
% 
% %% 
% 
% set(gcf,'Color','w')

%% also make the plot with the surrogate data that keeps cov


thrs = [50 95 99.9] ; 

tiledlayout(2,3)

for idx = 1:3
    nexttile()

    scatter(tv(spike_conn.subset1.long),tv(prctile(surrA.long,thrs(idx),3)),...
        'filled')
    
    axis square

    refline(1,0)

    if idx == 1
        ylabel({'randomized covariance' 'surrogate long event count'})
    end


    title(['surrogate percentile : ' num2str(thrs(idx))])

end


for idx = 1:3
    nexttile()

    scatter(tv(spike_conn.subset1.long),tv(prctile(surrB.long,thrs(idx),3)),...
        'filled')
    
    axis square

    refline(1,0)

    if idx == 1
        ylabel({'retained covariance' 'surrogate long event count'})
    end

    if idx == 2
        xlabel('empirical long event count')
    end

end

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[100 100 800 600])

%%

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_keepcov_scatters.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%% use the long edges for identifiability

scannames = {'REST1_LR' 'REST1_RL'} ; 

fmridat = struct() ; 
for idx = scannames
    fmridat.(idx{1}) = load_hcp_alldata(...
        [ DD.INTERM '/hcp352_nusregts_FIX2phys_schaefer200/' ],...
        'schaefer200-yeo17',...
        sublist.subset1, ...
        ['*' idx{1} '*']) ; 
end

%%
res.long = zeros(length(sublist.subset1),1) ; 
res.notlong = zeros(length(sublist.subset1),1) ; 

ident.long = zeros(length(sublist.subset1)) ; 
ident.notlong = zeros(length(sublist.subset1)) ; 

for idx = 1:length(sublist.subset1)

    disp(idx)

    sub = sublist.subset1{idx} ; 

    c1 = corr(fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes)) ; 
    c2 = corr(fmridat.REST1_LR(idx).ts(:,1:finfo.nnodes)) ; 

    res.long(idx) = corr(c1(triu(sigmoreA.long,1)),c2(triu(sigmoreA.long,1)),'type','s') ; 
    res.notlong(idx) = corr(c1(triu(~sigmoreA.long,1)),c2(triu(~sigmoreA.long,1)),'type','s') ; 

    % res.long(idx) = IPN_ccc([ c1(triu(sigmoreA.long,1)) c2(triu(sigmoreA.long,1)) ]) ; 
    % res.notlong(idx) = IPN_ccc([ c1(triu(~sigmoreA.long,1)) c2(triu(~sigmoreA.long,1)) ]) ; 

end

for idx = 1:length(sublist.subset1)
    disp(idx)
    c1 = corr(fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes)) ; 

    for jdx = 1:length(sublist.subset1)
        c2 = corr(fmridat.REST1_LR(jdx).ts(:,1:finfo.nnodes)) ; 

        ident.long(idx,jdx) = corr(c1(triu(sigmoreA.long,1)),c2(triu(sigmoreA.long,1)),'type','s') ; 
        ident.notlong(idx,jdx) = corr(c1(triu(~sigmoreA.long,1)),c2(triu(~sigmoreA.long,1)),'type','s') ; 
    end
end

% [~,mi1] = max(ident.long,[],2) ; 
% [~,mi2] = max(ident.notlong,[],2) ; 

[~,ss1] = sort(ident.long,2,'descend') ; 
identpos1 = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss,1)) ; 

[~,ss2] = sort(ident.notlong,2,'descend') ; 
identpos2 = arrayfun(@(i_) find(ss2(i_,:)==i_),1:size(ss,1)) ; 

dat = [ res.long,res.notlong ] ; 
[~,~,tt_ci,tt_stats] = ttest(dat(:,1),dat(:,2)) ; 

nperms = 1e4 ; 
permt = nan(nperms,1) ; 
rng(42)
for idx = 1:nperms
    dd = dat .* 1; 
    rr = randperm(176) ; 
    dd(rr(1:88),:) = fliplr(dd(rr(1:88),:)) ; 
    [~,~,~,ss] = ttest(dd(:,1),dd(:,2)) ; 
    permt(idx) = ss.tstat ; 
end

[rr_p,~,rr_stats] = ranksum(identpos1,identpos2) ;

%% 

tiledlayout(1,3)

nexttile()

imsc_grid_comm(sigmoreA.long,parc.ca(1:200),1,[0.5 0.5 0.5 ],[0.5 0.5 0.5 ],[])
axis square
colormap([ [0.8 0.8 0.8 ] ; CM(30,:) ])
xticks([])
yticks([])

ne = finfo.nnodes * (finfo.nnodes-1) / 2 ; 
legend({ ...
    ['sig. long (' num2str(round(sum(tv(sigmoreA.long))./ne*100,0)) '% edges) '] ...
    ['non-sig. long (' num2str(round(sum(tv(~sigmoreA.long))./ne*100,0)) '% edges)']}, ...
    'Location','southoutside')

nexttile()

h = histogram(res.long) ;
h.EdgeAlpha = 0 ; 
h.FaceColor = CM(30,:) ; 
h.FaceAlpha = 0.5 ; 
xline(mean(res.long),'Color',CM(30,:),'LineWidth',3)

hold on
h = histogram(res.notlong) ;
h.EdgeAlpha = 0 ; 
h.FaceColor = [0.8 0.8 0.8 ] ; 
h.FaceAlpha = 0.5 ; 
xline(mean(res.notlong),'Color',[.8 .8 .8],'LineWidth',3)

xlabel('test-retest correlation')
ylabel('count')

hold off

axis square

nexttile()

h,hh = fcn_boxpts([identpos1 identpos2]',...
    [ ones(176,1) ; ones(176,1)+1 ],[CM(30,:) ; [.8 .8 .8] ],...
    1,{'sig. long' 'non-sig. long'})
set(gca, 'YDir','reverse')

axis square

ylabel('identifiability rank')

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[100 100 800 600])

%%

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_vs_nonlongest_stats.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)


