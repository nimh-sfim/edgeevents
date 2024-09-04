%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%% 

thrs = 0:0.25:4 ; 
nthrs = length(thrs) ; 

res_nodewise = nan(NSUBS,finfo.nnodes,nthrs) ; 
res_glob = nan(NSUBS,nthrs) ; 

%%

nnodes = 200 ; 

for tdx = 1:nthrs

    disp([ num2str(tdx) '-' num2str(nthrs) ])
    tstart = tic ; 

    thr = thrs(tdx) ;

    tmpres_nw = nan(NSUBS,finfo.nnodes) ; 
    tmpres_glob = nan(NSUBS,1) ; 

    for idx = 1:NSUBS
        
        ts = datStr(idx).ts(:,1:nnodes) ; 
        ets = get_ets(ts) ; 
        fc = corr(ts) ; 

        spkc = count_spks(ets,thr) ; 

        tmpspk = mksq(spkc) ;  

        tmpres_nw(idx,:) = arrayfun(...
            @(i_)corr(fc(:,i_),tmpspk(:,i_),Type="Spearman"),...
            1:nnodes) ; 

        tmpres_glob(idx) = corr(tv(fc),tv(tmpspk),'Type','Spearman') ; 

    end

    res_nodewise(:,:,tdx) = tmpres_nw ; 
    res_glob(:,tdx) = tmpres_glob ; 

    disp(toc(tstart))

end

%% save these results

filename = [ DD.PROC '/spk_thr_dat_' , OUTSTR , '.mat'] ;
save(filename,'res_*','thrs')
 
%% 

TL = tiledlayout(1,2)

nexttile()

[~,mi] = max(res_glob,[],2) ; 
glob_thrs = thrs(mi) ; 

pickthr = mode(glob_thrs) ; % pickthr == 2

h = histogram(glob_thrs(:),BinEdges=thrs-0.125,EdgeAlpha=0.2) ; 
xlim([min(thrs) max(thrs)])
axis square
title({'global fc vs. spike count,' 'peak similarity values'})
xlabel('threshold values')
ylabel('count (subjects)')

text(0.05,0.05,[ 'mode: ' num2str(pickthr) ],'Units','normalized')

nexttile()

[~,mi] = max(res_nodewise,[],3) ; 
nodewise_thrs = thrs(mi) ; 

histogram(nodewise_thrs(:),BinEdges=thrs-0.125,EdgeAlpha=0.2)
xlim([min(thrs) max(thrs)])
axis square
title({'node-wise fc vs. spike count,' 'peak similarity values'})
xlabel('threshold values')
ylabel('count (nodes)')

pickthr = mode(nodewise_thrs(:)) ; % pickthr == 2

text(0.05,0.05,[ 'mode: ' num2str(pickthr) ],'Units','normalized')

% parc_plot(surfss,annotm,'schaefer200-yeo17',mean(nodewise_thrs),...
%         'cmap',parula(100),...
%         'viewcMap',1,'newFig',0,'viewStr','all')

% make a plot of the similarities at thr==2

out_figdir = [ './reports/figures/supp/' ] ; 
mkdir(out_figdir)
filename = [out_figdir '/spk_thr_sel_histograms.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%%

for idx =  [2 2.25]

    TL = tiledlayout(3,2)
    
    nt1 = nexttile([2 2])
    
    rrange = [0.6 0.9] ; 
    
    p = parc_plot(surfss,annotm,'schaefer200-yeo17',mean(res_nodewise(:,:,thrs==idx)),...
            'cmap',parula(100),'valRange',rrange,...
            'viewcMap',0,'newFig',0,'viewStr','all',...
            'parenth',TL) ;
    p.Layout = nt1.Layout ; 
    
    nt2 = nexttile(TL,[1 2])
    
    hh = imagesc(0:0.1:1,rrange) 
    cb = colorbar()
    hh.Visible = 'off' ;
    hh.Parent.Visible = 'off' ; 
    cb.Location = "north" ; 
    cm = colormap() ; 
    colormap(nt2,cm(2:100,:))
    cb.Label.String = 'fc vs. spike count, avg. similarity' ;
    
    %% write out the picture
    
    set(gcf,'Position',[100 100 1000 1000])
    set(gcf,'Color','w')

    out_figdir = [ './reports/figures/supp/' ]
    mkdir(out_figdir)
    filename = [out_figdir '/spk_thr_sel_map_thr' , num2str(idx) , '.png' ] ; 
    print(filename,'-dpng')
    close(gcf)

end