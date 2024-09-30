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
filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
load(filename)

%% boots

nboot = 1000 ; 

surrlong = zeros(finfo.nnodes,finfo.nnodes,nboot) ; 

nsubs = length(sublist.all) ; 
meansz = length(sublist.subset1) ; 

nperms = length(simmat.keepcov) ; 
for idx = 1:nboot
    
    disp(idx)

    % pickout the surr data to average
    bootinds = randsample(1:nperms,meansz,true) ; 

    tmp = zeros(finfo.nnodes) ; 
    for bdx = bootinds 

        ll = simmat.randcov{bdx}(3,:) ; 
        tmp = tmp + mksq(ll) ; 
    end
    
    surrlong(:,:,idx) = tmp./meansz ; 
    % surrlong(:,:,idx) = (tmp>prctile(tv(tmp),95)) ; 


end

%% sig test

% lowsim = prctile(surrlong,1,3) ; 
highsim = prctile(surrlong,99,3) ; 

sig_longedges = spike_conn.subset1.long>highsim ; 

%% 

dat = tv(spike_conn.subset1.long) ; 

TL = tiledlayout(1,2)

nt = nexttile()

longedge_thr = min(spike_conn.subset1.long(sig_longedges)) ; % prctile(dat,95) ; 
ll = linspace(0,longedge_thr,20) ; 
ss = ll(2) ; 

hh = histogram(dat,'BinEdges', ll ,FaceColor=[0.8 0.8 0.8],...
    Normalization='count')
hold on
histogram(dat(dat>=longedge_thr),'FaceColor',[0.2 0.2 0.2],...
    'BinEdges',ll(end):ss:max(dat),...
    Normalization='count')

ylabel('count')
xlabel('avg. num. long')

xx = xticks() ; 
xticklabels(cellstr(num2str(round(xx.*finfo.TR,2)')))

line([ longedge_thr longedge_thr ], [0 1750] ,'Color','red')
text(longedge_thr,2000, ...
    [ { 'surrogate' ['thr: ' num2str(round(longedge_thr*finfo.TR,2)) ] } ],...
    'HorizontalAlignment','center')

axis square

% aa = axes('Parent',gcf,'Position',[.48 .5 .4 .38],'Units','normalized')
% ss = scatter((dat1),(dat2),10,tld>3,"filled") ; 
% cc = inferno(10) ; 
% colormap(nt,[0.8 0.8 0.8 ; cc(2,:) ])
% rl = refline(1,0) ; 
% rl.LineWidth = 2 ;
% colormap(inferno(2))

nt1 = nexttile() 

longest_edges = dat >= longedge_thr ; 

bb = get_blocky(mksq(longest_edges),parc.ca(1:200))
imsc_grid_comm(bb,1:17,...
    1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],parc.names(1:17))
set(gca,'TickLength',[ 0 0])
colorbar
xticks([]) 

xlabel('btwn system density')

axis square

colormap(nt1,[ 1 1 1 ; flipud(inferno(100))])

set(gcf,'Position',[100 100 600 400])

%% 

set(gcf,'Color','w')

out_figdir = [ './reports/figures/figC/' ]
mkdir(out_figdir)
filename = [out_figdir '/longest_spike_hist.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)



%%

% nt2 = nexttile() 
% 
% pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', sum(mksq(longest_edges)) ,...
%     'valRange',[0 10],...
%     'cmap',[ 1 1 1 ; flipud(inferno(100))], ...
%     'viewcMap',0,'newFig',0,'viewStr','all',...
%     'parenth',TL)
% pp.Layout = nt2.Layout ; 


%%

addpath('/Users/faskowitzji/Documents/MATLAB/spm12/')
addpath('/Users/faskowitzji/Documents/MATLAB/conn/')

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

rois.sph_r = normalize(sum(mksq(longest_edges)), 'range',[1,5]) ; 

%%

viz_conn_glassbrain(mksq(longest_edges),[ 0.5 .5 .5],rois)

% 
% set(gcf,'Position',[100 100 600 400])
% 
% %% 
% 
% set(gcf,'Color','w')

