%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 
allspike_conn_indiv = struct() ; 
spkderiv = struct() ; 

for sdx = subsets

    allspike_conn_indiv.(sdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,length(sublist.(sdx{1}))) ; 
    spkderiv.(sdx{1}).mean = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spkderiv.(sdx{1}).burst1 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spkderiv.(sdx{1}).burst2 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    spkderiv.(sdx{1}).mean = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        spkmat = readdat.spike_len_mat  ; 
    
        [cc,dd] = count_spks(spkmat>0) ; 

        allspike_conn_indiv.(sdx{1})(:,:,idx) = mksq(cc) ; 

        % [a,~,c] = get_ctimes(spkmat>0) ; 
        % a = mksq(a) ; 
        % % b = mksq(b) ; 
        % c = mksq(c) ; 
        % 
        % iet_mats.(sdx{1}).mean(:,:,idx) = a(1:finfo.nnodes,1:finfo.nnodes) ; 
        % % iet_mats.(sdx{1}).max(:,:,idx) = b(1:finfo.nnodes,1:finfo.nnodes) ; 
        % iet_mats.(sdx{1}).std(:,:,idx) = c(1:finfo.nnodes,1:finfo.nnodes) ; 

        % a different mean calc
        % 
        % [a,~,c] = get_ctimes(dd) ;
        % a = mksq(a) ; 
        % % b = mksq(b) ; 
        % c = mksq(c) ; 

        % [burstVec,memVec,meanIET,maxIET] = spk_burst_mem_2(inSpk) 
        [bb1,bb2,mm,ee] = spk_burst_mem(dd) ; 

        bb1 = mksq(bb1) ; bb2 = mksq(bb2) ; mm = mksq(mm) ; ee = mksq(ee) ; 

        spkderiv.(sdx{1}).burst1(:,:,idx) = bb1(1:finfo.nnodes,1:finfo.nnodes) ; 
        spkderiv.(sdx{1}).burst2(:,:,idx) = bb2(1:finfo.nnodes,1:finfo.nnodes) ; 
        spkderiv.(sdx{1}).mem(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ; 
        spkderiv.(sdx{1}).mean(:,:,idx) = ee(1:finfo.nnodes,1:finfo.nnodes) ; 


    end

end
 
%%

filename = [ DD.PROC '/allspk_conn_indiv_' OUTSTR '.mat' ] ; 
save(filename,'allspike_conn_indiv','spkderiv','-v7.3')

%% and get the spike lengths 

subsets = {'subset1' 'subset2'} ; 
spike_len_mats = struct() ; 

for sdx = subsets

    spike_len_mats.(sdx{1}) = zeros(finfo.nnodes,finfo.nnodes,length(sublist.(sdx{1}))) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
       
        mm = mksq(cellfun(@mean,readdat.spike_len_cell )) ; 
        spike_len_mats.(sdx{1})(:,:,idx) = mm(1:finfo.nnodes,1:finfo.nnodes) ;
    end

end

%% save it

filename = [ DD.PROC '/allspk_len_indiv_' OUTSTR '.mat' ] ; 
save(filename,'spike_len_mats','-v7.3')

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

%% how many spikes per subject

eventtot = zeros(length(sublist.subset1),1) ; 
for idx = 1:length(sublist.subset1)

    eventtot(idx,1) = sum(tv(allspike_conn_indiv.subset1(1:200,1:200,idx))) ; 
    eventtot(idx,2) = sum(tv(allspike_conn_indiv.subset2(1:200,1:200,idx))) ; 
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

%% compare count vs fc

meanc = mean(allspike_conn_indiv.subset1,3) ; 
meanc = meanc(1:200,1:200) ; 

tiledlayout(1,3)

nexttile()
imsc_grid_comm(meanc,parc.ca(1:finfo.nnodes),0,[0.8 0.8 0.8],0.5)
colormap(flipud(blues()))
cb = colorbar ;
cb.Label.String = "count" ; 
axis square
xticks('')
yticks('')
title('mean count matrix')

nt = nexttile()
imsc_grid_comm(meanfc,parc.ca(1:finfo.nnodes),0,[0 0 0],0.5)
colormap(nt,rdbu())
cb = colorbar ;
cb.Label.String = "correlation mag." ; 
axis square
xticks('')
yticks('')
clim([-0.8 0.8])
title('mean correlation matrix')

nt = nexttile()
h = histscatter(tv(meanfc),tv(meanc),50) ;
% h.MarkerFaceColor = [0.5 0.2 0.5] ;
% h.MarkerFaceAlpha = 0.25 ; 
% h.MarkerEdgeAlpha = 0 ;
colormap(nt,purples())
cb = colorbar() ; 
cb.Label.String = 'bin count'
ylabel('count')
xlabel('correlation')
axis square
grid minor
title('mean corr. vs count')

%%

set(gcf,'Position',[100 100 800 200])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/corr_vs_count.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%

% compare count vs mean, per subject

countVmean1 = zeros(length(sublist.subset1),1) ; 
countVmean2 = zeros(length(sublist.subset1),1) ; 

for idx = 1:length(sublist.subset1)
    disp(idx)
    %mm = mksq(cellfun(@length,spike_lengths.subset1{idx})) ; 
    mm = allspike_conn_indiv.subset1(:,:,idx) ; 
    mm(isnan(mm)) = 0 ; 
    ff = fc_conn_indiv.subset1(:,:,idx) ; 
    countVmean1(idx) = corr(tv(mm),tv(ff),'type','s') ; 

    %mm = mksq(cellfun(@length,spike_lengths.subset2{idx})) ; 
    mm = allspike_conn_indiv.subset2(:,:,idx) ; 
    mm(isnan(mm)) = 0 ; 
    ff = fc_conn_indiv.subset2(:,:,idx) ; 
    countVmean2(idx) = corr(tv(mm),tv(ff),'type','s') ; 

end

%%

histogram(countVmean1,'BinWidth',.02)
hold on 
histogram(countVmean2,'BinWidth',.02)
hold off
xlim([0.5 1])
title('count vs. correlation, subject-level')
xlabel('rank correlation')
ylabel('num. subjects')

legend({'subset 1' 'subset 2'},'Location','northwest')

mean(countVmean1) 
std(countVmean1)
mean(countVmean2)
std(countVmean2)

% ans =
% 
%     0.8435
% 
% 
% ans =
% 
%     0.0429
% 
% 
% ans =
% 
%     0.8461
% 
% 
% ans =
% 
%     0.0461

%%

set(gcf,'Position',[100 100 300 300])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/corr_vs_count_indivhist.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%


out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/burst2_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%% and system maps

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
cmap = purples(25) ; 

fcn_boxpts(mean(countVcorr_nodewise)',...
    remap_labs,repmat(cmap(1,:),17,1),...
    0,parc.names(g1sort))
%([0 8.5])

set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

ylabel('mean coupling')
set(gca,"TickLabelInterpreter",'none')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_corrNcount_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% compare corrVcount across cortex

countVcorr_nodewise = zeros(length(sublist.subset1),finfo.nnodes) ; 
%countVcorr_nodewise2 = zeros(length(sublist.subset1),finfo.nnodes) ; 

for idx = 1:length(sublist.subset1)
    disp(idx)
    % mm = mksq(cellfun(@length,spike_lengths.subset1{idx})) ; 
    % mm(isnan(mm)) = 0 ; 
    mm = allspike_conn_indiv.subset1(:,:,idx) ; 
    ff = fc_conn_indiv.subset1(:,:,idx) ; 

    mm = mm(1:200,1:200); 
    ff = ff(1:200,1:200) ; 

    countVcorr_nodewise(idx,:) = arrayfun(@(i_) corr(mm(:,i_),ff(:,i_),'type','s'),1:finfo.nnodes) ;
    %countVcorr_nodewise2(idx,:) = arrayfun(@(i_) IPN_ccc([mm(:,i_) ff(:,i_)]),1:finfo.nnodes) ;

end

%%

parc_plot_wcolorbar(mean(countVcorr_nodewise),surfss,annotm,...
    [prctile(mean(countVcorr_nodewise),0) prctile(mean(countVcorr_nodewise),100)],flipud(purples(50)),[100 100 600 1000])

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
cmap = purples(25) ; 
%cmap = cmap(1:17,:) ; 


fcn_boxpts(mean(countVcorr_nodewise)',...
    remap_labs,repmat(cmap(1,:),17,1),...
    0,parc.names(g1sort))
%([0 8.5])

set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

ylabel('mean coupling')
set(gca,"TickLabelInterpreter",'none')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_corrNcount_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

np = 10000 ; 
pres = nan(np,1) ; 
for idx = 1:np 
    [~,tt,~] = anova1(mean(countVcorr_nodewise),parc.ca(randperm(200)),'off') ; 
    pres(idx) = tt{2,5} ; 
end
% emp
[~,ttt,sss] = anova1(mean(countVcorr_nodewise),parc.ca(1:200),'off') ; 
% ttt =
% 
%   4×6 cell array
% 
%     {'Source'}    {'SS'    }    {'df' }    {'MS'      }    {'F'       }    {'Prob>F'    }
%     {'Groups'}    {[0.1307]}    {[ 16]}    {[  0.0082]}    {[  7.6739]}    {[1.0791e-13]}
%     {'Error' }    {[0.1948]}    {[183]}    {[  0.0011]}    {0×0 double}    {0×0 double  }
%     {'Total' }    {[0.3255]}    {[199]}    {0×0 double}    {0×0 double}    {0×0 double  }

(sum(pres>=ttt{2,5})+1) / (np+1)

% tukey HSD

mm = multcompare(sss,'Alpha',0.001,'Display','off') ; 
%mmm = mksq(mm(:,6)) ; 
trilm = logical(tril(ones(17),-1)) ;
mmm = zeros(17) ; 
mmm(trilm) = mm(:,6) ;
mmm = mmm + mmm' ; 

mksq(fdr_bh(tv(mmm)))

%% finally, spike lengths

cmgreens = flipud(greens(100)) ; 

meanlen = mean(spike_len_mats.subset1,3,"omitmissing").*finfo.TR ; 
meanlen(1:size(meanlen)+1:end) = nan ;

parc_plot_wcolorbar(mean(meanlen,'omitmissing'),surfss,annotm,...
    [min(mean(meanlen,"omitmissing")) max(mean(meanlen,'omitmissing'))],cmgreens,[100 100 600 1000])

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/cortex_meanlength.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% compare count vs length

tiledlayout(1,3)

nt = nexttile()
imsc_grid_comm(meanlen,parc.ca(1:finfo.nnodes),0,[0 0 0],0.5)
colormap(nt,cmgreens)
cb = colorbar ;
cb.Label.String = "sec." ; 
axis square
xticks('')
yticks('')
title('mean event durration matrix')
bootci

nt = nexttile()
imsc_grid_comm(meanc,parc.ca(1:finfo.nnodes),0,[0.8 0.8 0.8],0.5)
colormap(nt,flipud(blues()))
cb = colorbar ;
cb.Label.String = "count" ; 
axis square
xticks('')
yticks('')
title('mean count matrix')

nt = nexttile()
h = histscatter(tv(meanc),tv(meanlen),50) ;
% h.MarkerFaceColor = [0.5 0.2 0.5] ;
% h.MarkerFaceAlpha = 0.25 ; 
% h.MarkerEdgeAlpha = 0 ;
colormap(nt,purples())
cb = colorbar() ; 
%clim([0 80])
cb.Label.String = 'bin count'
ylabel('durration')
xlabel('count')
axis square
grid minor
title('mean count vs durration')
corr(tv(meanc),tv(meanlen),'type','s')

%%

set(gcf,'Position',[100 100 800 200])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_vs_count.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% compare mean len (map) with other explanatory maps

ll = load('./data/interim/ts_variability.mat')

tsnrvec = (mean(ll.tsnrdat.REST1_LR.map,2) + mean(ll.tsnrdat.REST1_RL.map,2)) ./2 ; 
acvec = (mean(ll.acmap.REST1_LR.map,2) + mean(ll.acmap.REST1_RL.map,2)) ./2 ; 
dvarsvec = (mean(ll.dvarsdat.REST1_LR.map,2) + mean(ll.dvarsdat.REST1_RL.map,2)) ./2 ; 
cbvvec = squeeze(niftiread('./data/external/schaefer200_cbv-32k.pscalar.nii')) ; 
% cbfvec = squeeze(niftiread('./data/external/schaefer200_cbf-32k.pscalar.nii')) ; 

tiledlayout(2,2)

nexttile()
h = scatter_w_rho(mean(meanlen,2,'omitmissing'),acvec,'filled')
h.MarkerFaceColor = cmgreens(50,:) ; 
h.MarkerFaceAlpha = 0.5 ; 
h.MarkerEdgeColor = cmgreens(60,:) ; 

xlabel('node mean len. (sec)')
ylabel('autocor. (sec)')

h = refline(1,0)
h.Color = [0.8 0.8 0.8] ; 

axis square

% bootci(5000,{ @(a_,b_) corr(a_(:),b_(:),'type','s'), mean(meanlen,2,'omitmissing'), acvec},'type','per')
% 
% ans =
% 
%     0.7385
%     0.8565

nexttile()
h = scatter_w_rho(mean(meanlen,2,'omitmissing'),dvarsvec,'filled')
h.MarkerFaceColor = cmgreens(50,:) ; 
h.MarkerFaceAlpha = 0.5 ; 
h.MarkerEdgeColor = cmgreens(60,:) ; 

xlabel('node mean len. (sec)')
ylabel('dvars (a.u.)')

axis square

% bootci(5000,{ @(a_,b_) corr(a_(:),b_(:),'type','s'), mean(meanlen,2,'omitmissing'), dvarsvec},'type','per')
% ans =
% 
%    -0.4699
%    -0.2125

nexttile()
h = scatter_w_rho(mean(meanlen,2,'omitmissing'),tsnrvec,'filled')
h.MarkerFaceColor = cmgreens(50,:) ; 
h.MarkerFaceAlpha = 0.5 ; 
h.MarkerEdgeColor = cmgreens(60,:) ; 

xlabel('node mean len. (sec)')
ylabel('tsnr (a.u.)')

axis square

% bootci(5000,{ @(a_,b_) corr(a_(:),b_(:),'type','s'), mean(meanlen,2,'omitmissing'), tsnrvec},'type','per')
% ans =
%    -0.3475
%    -0.0736

nexttile()
h = scatter_w_rho(mean(meanlen,2,'omitmissing'),cbvvec,'filled')
h.MarkerFaceColor = cmgreens(50,:) ; 
h.MarkerFaceAlpha = 0.5 ; 
h.MarkerEdgeColor = cmgreens(60,:) ; 

xlabel('node mean len. (sec)')
ylabel('cerebral blood vol. (a.u.)')

axis square

% bootci(5000,{ @(a_,b_) corr(a_(:),b_(:),'type','s'), mean(meanlen,2,'omitmissing'), cbvvec},'type','per')
% ans =
% 
%   2×1 single column vector
% 
%    -0.1691
%     0.1303
set(gcf,'Position',[100 100 500 400])
set(gcf,'Color','w')

%% make system map

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

fcn_boxpts(mean(meanlen,2,'omitmissing'),...
    remap_labs,repmat(cmgreens(90,:),17,1),...
    0,parc.names(g1sort))
%([0 8.5])

set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

ylabel('node mean dur.')
set(gca,"TickLabelInterpreter",'none')

%%

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/meanlen_sys.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% and anova on that map

np = 10000 ; 
pres = nan(np,1) ; 
for idx = 1:np 
    [~,tt,~] = anova1(mean(meanlen,2,'omitmissing'),parc.ca(randperm(200)),'off') ; 
    pres(idx) = tt{2,5} ; 
end
% emp
[~,ttt,~] = anova1(mean(meanlen,2,'omitmissing'),parc.ca(1:200),'off') ; 
% ttt =
% 
%   4×6 cell array
% 
%     {'Source'}    {'SS'    }    {'df' }    {'MS'      }    {'F'       }    {'Prob>F'    }
%     {'Groups'}    {[2.0009]}    {[ 16]}    {[  0.1251]}    {[ 11.6768]}    {[1.4333e-20]}
%     {'Error' }    {[1.9599]}    {[183]}    {[  0.0107]}    {0×0 double}    {0×0 double  }
%     {'Total' }    {[3.9608]}    {[199]}    {0×0 double}    {0×0 double}    {0×0 double  }
(sum(pres>=ttt{2,5})+1) / (np+1)
% ans =
% 
%    9.9990e-05

%% get those highest areas

[~,si] = sort(mean(meanlen,2,'omitmissing'),'desc') ; 
annotm('schaefer200-yeo17').combo_names(si)

cortexplot(~cellfun(@isempty,regexpi(annotm('schaefer200-yeo17').combo_names,annotm('schaefer200-yeo17').combo_names(si(3)))))

%%

% lll = iet_mats.subset1.mean ; 
% lll(lll==finfo.ntp) = nan ; 
% 
% ietmean = mean(lll,3,'omitmissing') ; 
% 
% xdat = tv(meanc) ; 
% ydat = tv(ietmean) ; 
% 
% xbins = prctile(xdat,0:1:100) ; % just use simple percentiles, all bins have
%                                 % same numeber of edges
% [eVar,xbinmid,binmean] = norm_bin_model(xdat,xbins,ydat) ; 

%% TODO count vs corr figure??




% %% TODO -- will this be added to paper?
% 
% 
% %% let's look at the inter-event times
% 
% tmp = iet_mats.subset1.mean ; 
% tmp(tmp>=finfo.ntp-(floor(finfo.ntp*0.01))) = nan ; 
% ietmat = mean(tmp,3,'omitmissing') ; 
% 
% tmp = iet_mats.subset1.std ; 
% tmp(tmp>=finfo.ntp-(floor(finfo.ntp*0.01))) = nan ; 
% ietstdmat = mean(tmp,3,'omitmissing') ; 
% 
% 
% cmreds = flipud(reds()) ; 
% 
% tiledlayout(1,3) 
% 
% nt = nexttile() ; 
% 
% imsc_grid_comm(ietmat,parc.ca(1:finfo.nnodes),0,[0 0 0],0.5)
% colormap(nt,cmreds)
% cb = colorbar ;
% cb.Label.String = "sec." ; 
% axis square
% xticks('')
% yticks('')
% title('mean inter-event durr. matrix')
% 
% nt = nexttile()
% h = histscatter(tv(meanc),tv(ietmat),50) ;
% % h.MarkerFaceColor = [0.5 0.2 0.5] ;
% % h.MarkerFaceAlpha = 0.25 ; 
% % h.MarkerEdgeAlpha = 0 ;
% colormap(nt,purples())
% cb = colorbar() ; 
% %clim([0 80])
% cb.Label.String = 'bin count'
% ylabel('durration')
% xlabel('count')
% axis square
% grid minor
% title('mean count vs durration')
% corr(tv(meanc),tv(ietmat),'type','s')
% 
% nt = nexttile()
% h = histscatter(tv(meanlen),tv(ietmat),50) ;
% % h.MarkerFaceColor = [0.5 0.2 0.5] ;
% % h.MarkerFaceAlpha = 0.25 ; 
% % h.MarkerEdgeAlpha = 0 ;
% colormap(nt,purples())
% cb = colorbar() ; 
% %clim([0 80])
% cb.Label.String = 'bin count'
% ylabel('durration')
% xlabel('count')
% axis square
% grid minor
% title('mean count vs durration')
% corr(tv(meanlen),tv(ietmat),'type','s')
% 
% %%


%%

h = histscatter(tv(meanfc),tv(meanlen),50) ;
% h.MarkerFaceColor = [0.5 0.2 0.5] ;
% h.MarkerFaceAlpha = 0.25 ; 
% h.MarkerEdgeAlpha = 0 ;
colormap(purples())
cb = colorbar() ; 
%clim([0 80])
cb.Label.String = 'bin count'
ylabel('durration')
xlabel('correlation')
axis square
grid minor
title('correlation vs durration')
corr(tv(meanfc),tv(meanlen),'type','s')

set(gcf,'Position',[100 100 400 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figS/' ]
mkdir(out_figdir)
filename = [out_figdir '/corrVdurr_mat.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%



