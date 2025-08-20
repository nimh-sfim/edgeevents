%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 
%thrs = 0:0.25:4 ; 
thrs = (SPK_THR - (5*0.25)):0.25:(SPK_THR + (5*0.25)) ; 
nthrs = length(thrs) ; 

%% now read in spike lengths to make a histogram of lengths

%subsets = {'subset1' 'subset2'} ; 
subsets = {'subset1'} ; 
spkdurr = struct() ;

for sdx = subsets

    spkdurr.(sdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,nthrs) ; 

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) ' - ' sdx{1}])
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        tmpts = datStr(sind).ts ; 
        tmpets = get_ets(tmpts) ; 
        
        for tdx = 1:nthrs

            cthr = thrs(tdx) ; 

            abv_thr = tmpets > cthr ; 
   
            [~,dd] = spk_lenmat(abv_thr) ; 
            tmp = mksq(cellfun(@nanmean,dd)) ;
            tmp(isnan(tmp)) = 0 ;

            spkdurr.(sdx{1})(:,:,tdx) = spkdurr.(sdx{1})(:,:,tdx) + tmp ; 

        end

    end

    % after loop, average

    for tdx = 1:nthrs
        spkdurr.(sdx{1})(:,:,tdx) = spkdurr.(sdx{1})(:,:,tdx) ./ length(sublist.(sdx{1})) ; 
    end

end


filename = [ DD.PROC '/spk_durr_diffthrs_' , OUTSTR , '.mat'] ;
save(filename,'spkdurr','thrs','-v7.3')

%% quick sanity check analyses?

%load(filename)

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 

TL = tiledlayout(3,4)
TL.Title.String = 'group event duration matrices'

for idx = 1:nthrs
    nexttile()

    tmp = spkdurr.subset1(:,:,idx) .* finfo.TR ; 
    tmp = tmp(1:finfo.nnodes,1:finfo.nnodes) ; 

    imsc_grid_comm(tmp,...
        remap_labs,1,[0.2 0.2 0.2],[0.1 0.1 0.1 ],'')
    axis square
    %colorbar
    
    title(['thr: ' num2str(thrs(idx))])

    % cb = colorbar() ; 
    % cb.Label.String = 'average event count' ; 
    colormap("parula")
    % tmp = cb.TickLabels ; 
    % cb.TickLabels = { 'non sig.' tmp{2:end} } ; 
    cc = clim() ;  
    clim(cc) ;
    xticks([])
    yticks([])

    clim([1 4])
    
    if idx == nthrs
        cb = colorbar
        cb.Label.String = 'event duration (sec.)'
    end

end

colormap(cmgreens)

set(gcf,'Position',[100 100 800 500])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_thrs_allmats.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%
% res = nan(nthrs) ; 
% 
% for idx = 1:nthrs
%     tmp = spkdurr.subset1(:,:,idx) ; 
%     tmp = tmp(1:finfo.nnodes,1:finfo.nnodes) ; 
%     for jdx = 1:nthrs
%         tmp2 = spkdurr.subset1(:,:,jdx) ; 
%         tmp2 = tmp2(1:finfo.nnodes,1:finfo.nnodes) ; 
%         res(idx,jdx) = corr(tv(tmp),tv(tmp2)) ; 
%     end
% end

TL = tiledlayout(2,1)
TL.Title.String = 'comparison of event duration thresholds'
nexttile()

h = imsc_grid_comm(corr(cell2mat(arrayfun(@(i_) tv(spkdurr.subset1(:,:,i_)), 1:nthrs, 'UniformOutput', false)),'type','s') , ...
    1:nthrs,0.5,[1 1 1],[1 1 1],cellstr(string(thrs)))
axis square
%colormap(cmgreens)
cb = colorbar ;
cb.Label.String = 'rank correlation' ; 
xticks('')
ylabel('event thresholds')
h.AlphaData = ~eye(nthrs)
set(gca,'color',[0.5 0.5 0.5])
set(gcf,'color',[1 1 1])

nexttile()

h = imsc_grid_comm(lins_ccc(cell2mat(arrayfun(@(i_) tv(spkdurr.subset1(:,:,i_)), 1:nthrs, 'UniformOutput', false))) , ...
    1:nthrs,0.5,[1 1 1],[1 1 1],cellstr(string(thrs)))
axis square
%colormap(cmgreens)
cb = colorbar ;
cb.Label.String = "lin's concordance"
xticks('')
h.AlphaData = ~eye(nthrs)
set(gca,'color',[0.5 0.5 0.5])
set(gcf,'color',[1 1 1])

colormap("gray")

%%

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_thrs_mats.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

% %% make the cortex plots for all the different thresholds?
% 
% for idx = 1:11
% 
%     dd = spkdurr.subset1(:,:,idx) ; 
%     dd = dd(1:finfo.nnodes,1:finfo.nnodes) ; 
% 
%     % parc_plot_wcolorbar(mean(dd,2),surfss,annotm,...
%     %     [min(mean(dd,2)) max(mean(dd,2))],cmgreens,[100 100 600 1000])
% 
%     imagesc(dd)
%     title(idx)
% 
%     waitforbuttonpress
% 
%     % out_figdir = [ './reports/figures/figR/' ] ; 
%     % mkdir(out_figdir)
%     % orient(gcf,'landscape')
%     % filename = [out_figdir '/durr_' datnames{idx} '_cortex.pdf' ] ; 
%     % print(filename,'-dpdf')
%     % close(gcf)
% 
% end

%%

filename = [ DD.PROC '/spk_thr_dat_' , OUTSTR , '.mat'] ;
ll = load(filename)

%%

filename = [ DD.PROC '/allspk_conn_indiv_' OUTSTR '.mat' ] ; 
ss = load(filename,'allspike_conn_indiv','spkderiv') ; 


%%

% quick get the inds for subset 1
subset1pickout = arrayfun(@(i_)find(cellfun(@(x_)strcmp(x_,sublist.subset1(i_)),sublist.all)),1:176) ; 
tt = ll.thrs ; 

[~,mi] = max(ll.res_glob(subset1pickout,:),[],2) ;
umi = unique(mi) ; 
res = nan(finfo.nnodes,finfo.nnodes,length(umi)) ; 


tiledlayout(2,4)

for idx = 1:length(umi)

    nexttile()

    uu = umi(idx) ; 
    pickout = mi==uu ; 

    tmp = ss.allspike_conn_indiv.subset1(:,:,pickout) ; 
    res(:,:,idx) = mean(tmp(1:finfo.nnodes,1:finfo.nnodes,:),3) ; 

    imsc_grid_comm(res(:,:,idx),...
        remap_labs,[],[],[],'')
    xticks('') ; yticks('')
    axis square
    cb = colorbar ;
    clim([0 30])
        
    title([ 'thr: ' num2str(tt(uu)) ', #subs: ' num2str(sum(pickout)) ])
    %title([ titlestrings{idx} ', num. subs: ' num2str(sum(pickout)) ])

    if idx == length(umi)
        cb.Label.String = 'num. events'
    end

end
 
set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figR/' ]
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_spk_mats.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% and now the correlation

%ccc = lins_ccc(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false))) ; 
ccc = corr(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false)),'type','s') ; 

h = imsc_grid_comm(ccc , ...
    1:length(umi),0.5,[1 1 1],[1 1 1],cellstr(string(tt(umi))))
axis square
colormap("parula")
cb = colorbar ;
cb.Label.String = 'rank correlation' ; 
xticks('')
h.AlphaData = ~eye(length(umi))
set(gca,'color',[0.5 0.5 0.5])
set(gcf,'color',[1 1 1])

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

clim([0 1])

out_figdir = [ './reports/figures/figR/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_spk_mats_corr.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

ccc = lins_ccc(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false))) ; 
%ccc = corr(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false)),'type','s') ; 

h = imsc_grid_comm(ccc , ...
    1:length(umi),0.5,[1 1 1],[1 1 1],cellstr(string(tt(umi))))
axis square
colormap("parula")
cb = colorbar ;
cb.Label.String = "lin's concord." ; 
xticks('')
h.AlphaData = ~eye(length(umi))
set(gca,'color',[0.5 0.5 0.5])
set(gcf,'color',[1 1 1])

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figR/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_spk_mats_lins.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%%

%%

% quick get the inds for subset 1
subset1pickout = arrayfun(@(i_)find(cellfun(@(x_)strcmp(x_,sublist.subset1(i_)),sublist.all)),1:176) ; 
tt = ll.thrs ; 

[~,mi] = max(ll.res_glob(subset1pickout,:),[],2) ;
umi = unique(mi) ; 
mi = dummyvaroppo([ ismember(mi,umi(tt(unique(mi))<SPK_THR)) ismember(mi,umi(tt(unique(mi))==SPK_THR)) ismember(mi,umi(tt(unique(mi))>SPK_THR)) ]) ; 

umi = unique(mi) ; 
res = nan(finfo.nnodes,finfo.nnodes,length(umi)) ; 

titlestrings = { '<2.25' '=2.25' '>2.25'} ; 

tiledlayout(1,3)

for idx = 1:length(umi)

    nexttile()

    uu = umi(idx) ; 
    pickout = mi==uu ; 

    tmp = ss.allspike_conn_indiv.subset1(:,:,pickout) ; 
    res(:,:,idx) = mean(tmp(1:finfo.nnodes,1:finfo.nnodes,:),3) ; 

    imsc_grid_comm(res(:,:,idx),...
        remap_labs,[],[],[],'')
    xticks('') ; yticks('')
    axis square
    cb = colorbar ;
    clim([0 30])
        
    title([ titlestrings{idx} ', num. subs: ' num2str(sum(pickout)) ])

    if idx == length(umi)
        cb.Label.String = 'num. events'
    end

end
 
set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/figR/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_spk_mats2.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%% 
% 

% [~,mi] = max(ll.res_glob(subset1pickout,:),[],2) ;
% umi = unique(mi) ; 
% 
% niter = 1000 ; 
% resiter = nan(length(umi),length(umi),niter) ; 
% tmpres = nan(finfo.nnodes,finfo.nnodes,length(umi)) ; 
% 
% for iidx = 1:niter
% 
%     disp(iidx)
% 
%     randmi = mi(randperm(length(mi))) ; 
% 
%     for idx = 1:length(umi)
% 
%         uu = umi(idx) ; 
%         pickout = randmi==uu ; 
% 
%         tmp = ss.allspike_conn_indiv.subset1(:,:,pickout) ; 
%         tmpres(:,:,idx) = mean(tmp(1:finfo.nnodes,1:finfo.nnodes,:),3) ; 
%     end
% 
%     resiter(:,:,iidx) = ...
%         corr(cell2mat(arrayfun(@(i_) tv(tmpres(:,:,i_)), 1:length(umi), 'UniformOutput', false)),'type','s') ;
% 
%     % resiter(:,:,iidx) = ...
%     %     lins_ccc(cell2mat(arrayfun(@(i_) tv(tmpres(:,:,i_)), 1:length(umi), 'UniformOutput', false))) ;
% 
% end
% 
% ccc = corr(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false)),'type','s') ; 


%%

ccc = corr(cell2mat(arrayfun(@(i_) tv(res(:,:,i_)), 1:length(umi), 'UniformOutput', false)),'type','s') ; 

% ccc =
% 
%     1.0000    0.9750    0.9597
%     0.9750    1.0000    0.9718
%     0.9597    0.9718    1.0000


h = imsc_grid_comm( ccc , ...
    1:length(umi),0.5,[1 1 1],[1 1 1],cellstr(string(tt(umi))))
axis square
colormap("parula")
cb = colorbar ;
cb.Label.String = 'rank correlation' ; 
xticks('')
h.AlphaData = ~eye(length(umi))
set(gca,'color',[0.5 0.5 0.5])
set(gcf,'color',[1 1 1])

clim([0 1])

set(gcf,'Position',[100 100 1000 400])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figR/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/diff_durr_spk_mats_corr2.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%%

dd = std(meanlen,'omitmissing') ; 
parc_plot_wcolorbar(dd,surfss,annotm, ...
    [min(dd) max(dd)],cmgreens,[100 100 600 1000])


dd = mean(mean(ss.spkderiv.subset1.burst1,3,'omitmissing')) ; 
parc_plot_wcolorbar(dd,surfss,annotm, ...
    [min(dd) max(dd)],cmgreens,[100 100 600 1000])

%% figure our the range of duration mean relative to all durations 

% get all durs at 2.25 thr

%subsets = {'subset1' 'subset2'} ; 
subsets = {'subset1'} ; 
spkdurr2p25 = zeros(finfo.nnodes,finfo.nnodes,length(sublist.subset1)) ; 

for idx = 1:length(sublist.subset1)

    disp([ num2str(idx) ' - ' ])

    sind = find(cellfun(@(x_)strcmp(x_,sublist.subset1(idx)),sublist.all)) ; 

    tmpts = datStr(sind).ts ; 
    tmpets = get_ets(tmpts) ; 
    
    abv_thr = tmpets > SPK_THR ; 

    [~,dd] = spk_lenmat(abv_thr) ; 
    tmp = mksq(cellfun(@nanmean,dd)) ;
    tmp(isnan(tmp)) = 0 ;
    tmp = tmp(1:200,1:200) ; 
        
    spkdurr2p25(:,:,idx) =  tmp ; 

end

min(cell2mat(arrayfun(@(i_) mean(spkdurr2p25(:,:,i_)) , 1:176 , 'UniformOutput' , false)')).*finfo.TR ; 
max(cell2mat(arrayfun(@(i_) mean(spkdurr2p25(:,:,i_)) , 1:176 , 'UniformOutput' , false)')).*finfo.TR ; 

std(cell2mat(arrayfun(@(i_) mean(spkdurr2p25(:,:,i_),2) , 1:176 , 'UniformOutput' , false)')).*finfo.TR

