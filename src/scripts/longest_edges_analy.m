%% start
%% TODO ADAPT THIS FOR SIG TESTING ALL EDGES 

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ;

%% load the spike conn

filename = [ DD.PROC '/spk_conn_ex_avg_' OUTSTR '.mat' ] ; 
load(filename)

%% do some stat testing

% load up the nulls
filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
load(filename)

%% boots

nboot = 5000 ; 

nsubs = length(sublist.all) ; 
subsetsz = length(sublist.subset1) ; 

nperms = length(simmat.keepcov) ; 
rng(42)

surrA = struct() ;
surrB = struct() ;

lennames = {'longer' 'longest'} ; 

% and record for this boot
for sdx = 1:2
    surrA.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot); 
    surrB.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot) ; 
end

for idx = 1:nboot
    
    disp(idx)

    % pickout the surr data to average
    bootinds = randsample(1:nperms,subsetsz,true) ; 

    tmpA = zeros(2,finfo.nnodes*(finfo.nnodes-1)/2) ; 
    tmpB = zeros(2,finfo.nnodes*(finfo.nnodes-1)/2)  ; 

    for bdx = bootinds 

        aa = arrayfun(@(i_)simmat.nocov{bdx}(i_,:),4:5,'UniformOutput',false) ; 
        bb = arrayfun(@(i_)simmat.keepcov{bdx}(i_,:),4:5,'UniformOutput',false) ; 

        for sdx = 1:2
            tmpA(sdx,:) = tmpA(sdx,:) + aa{sdx}  ; 
            tmpB(sdx,:) = tmpB(sdx,:) + bb{sdx}  ; 
        end

    end
    
    if any(any(mksq(tmpB(1,:))<mksq(tmpB(2,:))))
        error('whats up here')
    end

    % make mean
    for sdx = 1:2
        tmpA(sdx,:) = tmpA(sdx,:) ./ subsetsz  ; 
        tmpB(sdx,:) = tmpB(sdx,:) ./ subsetsz  ; 
    end

    % and record for this boot
    for sdx = 1:2
        surrA.(lennames{sdx})(:,:,idx) = mksq(tmpA(sdx,:)) ; 
        surrB.(lennames{sdx})(:,:,idx) = mksq(tmpB(sdx,:)) ; 
    end

end 

%%

% sig = struct() ; 
% 
% for sdx = 1:2
% 
%     sig.A.more.(lennames{sdx}) = spike_conn_ex.subset1.(lennames{sdx}) >= ...
%         prctile(surrA.(lennames{sdx}),100,3) ; 
%     sig.B.more.(lennames{sdx}) = spike_conn_ex.subset1.(lennames{sdx}) >= ...
%         prctile(surrB.(lennames{sdx}),100,3) ; 
% 
%     sig.A.less.(lennames{sdx}) = spike_conn_ex.subset1.(lennames{sdx}) <= ...
%         prctile(surrA.(lennames{sdx}),0,3) ; 
%     sig.B.less.(lennames{sdx}) = spike_conn_ex.subset1.(lennames{sdx}) <= ...
%         prctile(surrB.(lennames{sdx}),0,3) ; 
% 
% end

%% 

sig = struct() ; 

for sdx = 1:2

    sig.A.more.(lennames{sdx}) = mksq(fdr_bh(tv((sum(surrA.(lennames{sdx})>=spike_conn_ex.subset1.(lennames{sdx}),3)+1) ./ (nboot+1)))) ; 
    sig.B.more.(lennames{sdx}) = mksq(fdr_bh(tv((sum(surrB.(lennames{sdx})>=spike_conn_ex.subset1.(lennames{sdx}),3)+1) ./ (nboot+1)))) ; 

    sig.A.less.(lennames{sdx}) = mksq(fdr_bh(tv((sum(surrA.(lennames{sdx})<=spike_conn_ex.subset1.(lennames{sdx}),3)+1) ./ (nboot+1)))) ; 
    sig.B.less.(lennames{sdx}) = mksq(fdr_bh(tv((sum(surrB.(lennames{sdx})<=spike_conn_ex.subset1.(lennames{sdx}),3)+1) ./ (nboot+1)))) ; 

end


%% save the sig info

filename = [ DD.PROC '/spk_longestedge_sigmask_' OUTSTR '.mat' ] ; 
save(filename,'sig','-v7.3')

%% plot eeet

CM = flipud(blues(100)) ; 

TL = tiledlayout(3,2) ;
set(gcf,'Position',[100 100 800 800])

nt = nexttile() ;

dat = spike_conn_ex.subset1.longest ; 

imsc_grid_comm(dat,remap_labs,1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],parc.names(g1sort))
axis square

cb = colorbar() ; 
cb.Label.String = 'average event count' ; 
colormap(nt,CM)
cc = clim() ;  
clim(cc) ;
xticks([])


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

nt = nexttile()

sigmask = sig.B.more.longest ; 
dat = spike_conn_ex.subset1.longest ; 

imsc_grid_comm(sigmask,remap_labs,1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],parc.names(g1sort))
axis square

% cb = colorbar() ; 
% cb.Label.String = 'average event count' ; 
colormap(nt,[ 1 1 1 ; CM(80,:)])
% tmp = cb.TickLabels ; 
% cb.TickLabels = { 'non sig.' tmp{2:end} } ; 
cc = clim() ;  
clim(cc) ;
xticks([])


nt2 = nexttile() ;

%longest_edges = dat >= longedge_thr ; 

bb = get_blocky(sigmask,remap_labs)
imsc_grid_comm(bb,1:17,...
    1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],[])
set(gca,'TickLength',[ 0 0])
cb = colorbar()
cb.Label.String = 'density of sig. edges' ; 
xticks([]) 
yticks([])

% xlabel('btwn system density')

axis square

colormap(nt2,[ 1 1 1 ; CM])

nexttile()

histogram(nonzeros(tv(sigmask.*dat)),...
    'Normalization','count','FaceColor',CM(80,:),'EdgeAlpha',0,'BinWidth',0.025)
hold on
histogram(nonzeros(tv(~sigmask.*dat)),...
    'Normalization','count','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0,'BinWidth',0.025)

%legend({'sig. edges' 'non-sig. edges'},'Location','northeastoutside')

xlabel('average event count')
ylabel('num. edges')

axis square

nexttile()

histogram(nonzeros(tv(sigmask.*meanfc)),...
    'Normalization','count','FaceColor',CM(80,:),'EdgeAlpha',0,'BinWidth',0.025)
xline(mean(nonzeros(tv(sigmask.*meanfc))),Color=CM(80,:))

hold on
histogram(nonzeros(tv(~sigmask.*meanfc)),...
    'Normalization','count','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0,'BinWidth',0.025)

xline(mean(nonzeros(tv(~sigmask.*meanfc))),Color=[0.5 0.5 0.5])


legend({'sig. edges' 'non-sig. edges'},'Location','northeastoutside')

xlabel('correlation (FC)')
ylabel('num. edges')

axis square



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
% % rois.sph_r = normalize(sum(sigmoreA.long), 'range',[1,3]) ; 
% % viz_conn_glassbrain(sigmoreA.long,[ 0.5 .5 .5],rois)
% % 
% out_figdir = [ './reports/figures/figC/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/longest_spike_glassbrain.png' ] ; 
% h('print',3,filename,'-nogui')
% close(gcf)

parc_plot_wcolorbar(sum(sigmask),surfss,annotm,...
    [0 max(sum(sigmask))],CM,[100 100 600 1000])

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longestsig_cortex.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)


% 
% set(gcf,'Position',[100 100 600 400])
% 
% %% 
% 
% set(gcf,'Color','w')

%% breakdown by system

cmap = CM(end-25:end-9,:) ; 

fcn_boxpts(sum(sigmask,2),...
    remap_labs,cmap,...
    0,parc.names(g1sort))
%([0 8.5])

%ylim([ 0 .1])
axis square
ylabel('degree')

%%

set(gca, 'TickLabelInterpreter', 'none')
set(gcf,'Position',[100 100 1300 400])

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_spike_by_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% some stats

% function [anovastuff,tbl,stats] = perm_mat_anova1(mat,labs,nperms)
%[aa,bb,cc] = perm_mat_anova1(sigmask,remap_labs,100) ; 
np = 10000 ; 
pres = nan(np,1) ; 
for idx = 1:np 
    [~,tt,~] = anova1(sum(sigmask),remap_labs(randperm(length(remap_labs))),'off') ; 
    pres(idx) = tt{2,5} ; 
end
% emp
[~,ttt,~] = anova1(sum(sigmask),remap_labs,'off') ; 
% ttt =
% 
%   4×6 cell array
% 
%     {'Source'}    {'SS'        }    {'df' }    {'MS'        }    {'F'       }    {'Prob>F'    }
%     {'Groups'}    {[1.6783e+05]}    {[ 16]}    {[1.0490e+04]}    {[ 21.7321]}    {[3.7892e-34]}
%     {'Error' }    {[8.8330e+04]}    {[183]}    {[  482.6760]}    {0×0 double}    {0×0 double  }
%     {'Total' }    {[2.5616e+05]}    {[199]}    {0×0 double  }    {0×0 double}    {0×0 double  }

(sum(pres>=ttt{2,5})+1)/(np+1)
% ans =
% 
%    9.9990e-05

%% also make the plot with the surrogate data that keeps cov


thrs = [50 95 99.9] ; 

tiledlayout(2,3)

for idx = 1:3
    nexttile()

    scatter(tv(spike_conn_ex.subset1.longest),tv(prctile(surrA.longest,thrs(idx),3)),...
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

    scatter(tv(spike_conn_ex.subset1.longest),tv(prctile(surrB.longest,thrs(idx),3)),...
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


%% load time series

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

sigmask = sig.B.more.longest; 

%ll = load()

for idx = 1:length(sublist.subset1)

    disp(idx)

    sub = sublist.subset1{idx} ; 

    c1 = corr(fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes)) ; 
    c2 = corr(fmridat.REST1_LR(idx).ts(:,1:finfo.nnodes)) ; 

    % res.long(idx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','s') ; 
    % res.notlong(idx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','s') ; 

    res.long(idx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
    %res.notlong(idx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 

    nperm=100 ; 
    tmpres = nan(nperm,1) ; 
    for jdx = 1:nperm
        % shuffmask = ~~randomizer_bin_und(sigmask,100) ; 
        shuffmask = randmio_und(sigmask,5) ; 
        tmpres(jdx) = IPN_ccc([ c1(triu(shuffmask,1)) c2(triu(shuffmask,1)) ]) ; 
    end
    res.notlong(idx) = mean(tmpres,'omitnan') ;
end

[~,tt_p,tt_ci,tt_stats] = ttest(res.long(:),res.notlong(:)) ; 

%%

rng(42)
fcn_boxpts([res.long(:) ; res.notlong(:)], ...
    [ones(176,1) ; ones(176,1).*2],[CM(80,:) ; 0.5 0.5 0.5 ],1,{'sig.','non-sig'})
ylabel("Lin's concordance")

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_testretest.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%% doesn't work
% for idx = 1:length(sublist.subset1)
%     disp(idx)
%     c1 = fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes) ; 
%     c1ets = get_ets(c1) ; 
% 
% 
% 
%     for jdx = 1:length(sublist.subset1)
%         c2 = corr(fmridat.REST1_LR(jdx).ts(:,1:finfo.nnodes)) ; 
% 
%         ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
%         ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
%         % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
%         % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
%     end
% end
% 
% % [~,mi1] = max(ident.long,[],2) ; 
% % [~,mi2] = max(ident.notlong,[],2) ; 
% %%
% 
% [~,ss1] = sort(ident.long,2,'descend') ; 
% identpos1 = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 
% 
% [~,ss2] = sort(ident.notlong,2,'descend') ; 
% identpos2 = arrayfun(@(i_) find(ss2(i_,:)==i_),1:size(ss1,1)) ; 

%%

% filename = [ DD.PROC '/spk_conn_ex_indiv_' OUTSTR '.mat' ] ; 
% ll1 = load(filename) ; 
% 
% filename = [ DD.PROC '/spk_conn_ex_indiv_' OUTSTR '_2.mat' ] ; 
% ll2 = load(filename) ; 

%%

% for idx = 1:length(sublist.subset1)
%     disp(idx)
%     c1 = corr(fmridat.REST1_RL(idx).ts(:,1:finfo.nnodes)) ; 
% 
%     for jdx = 1:length(sublist.subset1)
%         c2 = corr(fmridat.REST1_LR(jdx).ts(:,1:finfo.nnodes)) ; 
% 
%         ident.long(idx,jdx) = corr(c1(triu(sigmask,1)),c2(triu(sigmask,1)),'type','p') ; 
%         ident.notlong(idx,jdx) = corr(c1(triu(~sigmask,1)),c2(triu(~sigmask,1)),'type','p') ;
%         % ident.long(idx,jdx) = IPN_ccc([ c1(triu(sigmask,1)) c2(triu(sigmask,1)) ]) ; 
%         % ident.notlong(idx,jdx) = IPN_ccc([ c1(triu(~sigmask,1)) c2(triu(~sigmask,1)) ]) ; 
%     end
% end
% 
% [~,mi1] = max(ident.long,[],2) ; 
% [~,mi2] = max(ident.notlong,[],2) ; 


%%

[~,ss1] = sort(ident.long,2,'descend') ; 
identpos1 = arrayfun(@(i_) find(ss1(i_,:)==i_),1:size(ss1,1)) ; 

[~,ss2] = sort(ident.notlong,2,'descend') ; 
identpos2 = arrayfun(@(i_) find(ss2(i_,:)==i_),1:size(ss1,1)) ; 


%%

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

%[rr_p,~,rr_stats] = ranksum(identpos1,identpos2) ;

%% 

tiledlayout(1,2)

nexttile()

sigmask = sig.B.more.longest ; 

imsc_grid_comm(sigmask,parc.ca(1:200),1,[0.5 0.5 0.5 ],[0.5 0.5 0.5 ],[])
axis square
colormap([ [0.8 0.8 0.8 ] ; CM(30,:) ])
xticks([])
yticks([])

ne = finfo.nnodes * (finfo.nnodes-1) / 2 ; 
legend({ ...
    ['sig. long (' num2str(round(sum(tv(sigmask))./ne*100,0)) '% edges) '] ...
    ['non-sig. long (' num2str(round(sum(tv(~sigmask))./ne*100,0)) '% edges)']}, ...
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

xlabel('test-retest concordance')
ylabel('count')

hold off

axis square

% nexttile()
% 
% h,hh = fcn_boxpts([identpos1 identpos2]',...
%     [ ones(176,1) ; ones(176,1)+1 ],[CM(30,:) ; [.8 .8 .8] ],...
%     1,{'sig. long' 'non-sig. long'})
% set(gca, 'YDir','reverse')
% 
% axis square
% 
% ylabel('identifiability rank')
% 
% set(gcf,'Color','w')
% orient(gcf,'landscape')
% set(gcf,'Position',[100 100 800 600])

%%

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_vs_nonlongest_stats.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

% %% does knowing the fc annnnd the long spikes help out??
% 
% 
% filename = [ DD.PROC '/spk_conn_indiv_' OUTSTR '.mat' ] ; 
% load(filename)
% 
% filename = [ DD.PROC '/fc_conn_indiv_' OUTSTR '.mat' ] ; 
% load(filename)

%%




