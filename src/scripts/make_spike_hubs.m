%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

high_bin = 4 ; 
maxspk = 1200  ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

spklen_names = {'short' 'inter' 'long'} ; 

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

refparc = parc.ca(1:200) ; 

%% now read in spike lengths to make a histogram of lengths

spikestr = struct() ; 

mat_mst = @(x_) full(adjacency(minspantree(graph(x_)))) ; 

for sdx = subsets

    for hdx = 1:3
        spikestr.dd.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
        spikestr.ss.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
        spikestr.ec.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
        spikestr.cc.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
        spikestr.cl.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
        spikestr.bb.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 

    end

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        %origfc = corr(datStr(sind).ts(:,1:200)) ; 
        %mmask = threshold_proportional(corr(origfc),0.2) ; 
        %mmask = percthr_trimlow_und(origfc) ; 

        %% lets figure out if the hubs are the same! 

        for hdx = 1:3

            mat = mksq(mean(dd==hdx)) ; 

            % [mmask,ii] = my_omst(-(mat-max(mat,[],'all'))) ; 
            % mmask = threshold_proportional(mat,0.2) ; 
            mmask = percthr_trimlow_und(mat) ; 

            spikestr.dd.(sdx{1}).(spklen_names{hdx})(:,idx) = sum(mmask) ;  
            spikestr.cc.(sdx{1}).(spklen_names{hdx})(:,idx) = clustering_coef_bu(mmask) ; 
            
            spikestr.ss.(sdx{1}).(spklen_names{hdx})(:,idx) = sum(mat.*mmask) ;  
            spikestr.ec.(sdx{1}).(spklen_names{hdx})(:,idx) = eigenvector_centrality_und(mat.*mmask) ; 

            spikestr.cl.(sdx{1}).(spklen_names{hdx})(:,idx) = my_close_cent(mat.*mmask) ; 
            spikestr.bb.(sdx{1}).(spklen_names{hdx})(:,idx) = my_btwn_cent(mat.*mmask) ; 

        end

    end
end

%% save it

filename = [ DD.PROC '/spk_hubs_' OUTSTR '.mat' ] ; 
save(filename,'spikestr','-v7.3')


%%

figure
cortexplot(mean(spikestr.zz.subset2.(spklen_names{1}),2))

%% divide out by yeo sys

% g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
% 
% remap_nodes = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
% [datplot_ca,datplot_sinds] = sort(remap_nodes) ;  
% datplot_lab_sinds = remaplabs(1:17,g1sort,1:17) ; 
% cmap = get_nice_yeo_cmap('grad1') ; 

% lims = [50 150] ; 

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

% clf
% TL = tiledlayout(1,3,'TileIndexing','columnmajor') ; 
CM = flipud(gnbu(100)) ; 

for idx = 1:3

    % nt1 = nexttile(TL)

    % combo = [ mean(spikestr.dd.subset1.(spklen_names{idx}),2,'omitnan') ...
    %     mean(spikestr.cc.subset1.(spklen_names{idx}),2,'omitnan') ...
    %     mean(spikestr.ss.subset1.(spklen_names{idx}),2,'omitnan') ...
    %     mean(spikestr.ec.subset1.(spklen_names{idx}),2,'omitnan') ...
    %     mean(spikestr.cl.subset1.(spklen_names{idx}),2,'omitnan') ...
    %     mean(spikestr.bb.subset1.(spklen_names{idx}),2,'omitnan') ] ; 
    combo = zscore([ mean(normalize(spikestr.dd.subset1.(spklen_names{idx}),'range',[0 1]),2) ...
        mean(normalize(-spikestr.cc.subset1.(spklen_names{idx}),'range',[0 1]),2) ...
        mean(normalize(spikestr.bb.subset1.(spklen_names{idx}),'range',[0 1]),2) ]) ; 
    [~,dat] = pca(combo,'NumComponents',1) ; 


    %dat = mean(spikestr.cc.subset1.(spklen_names{idx}),2,'omitnan') ;
    lims = [min(dat) max(dat)] ; 

    parc_plot_wcolorbar(dat,surfss,annotm,...
        lims,CM,[100 100 600 1000])

    set(gcf,'Color','w')
    
    out_figdir = [ './reports/figures/figA/' ]
    mkdir(out_figdir)
    filename = [out_figdir '/spike_hub_' spklen_names{idx} '.png' ] ; 
    print(filename,'-dpng')
    close(gcf)

    % pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dat ,...
    %     'valRange',lims,...
    %     'cmap',CM, ...
    %     'viewcMap',0,'newFig',0,'viewStr','all',...
    %     'parenth',TL)
    % pp.Layout = nt1.Layout ; 
    % 
    % nt1.Title.String = spklen_names{idx} ; 
    % nt1.Title.Visible = 'on' ; 

    % nexttile(TL)
    % 
    % % fcn_boxpts(dat(datplot_sinds),...
    % %     datplot_ca,cmap(datplot_lab_sinds,:),...
    % %     0,parc.names(g1sort))
    % % 
    % % ylim(lims)
    % 
    % scatter(grad_cifti(1,:),grad_cifti(2,:),40,dat,'filled','MarkerEdgeColor','flat')
    % xlabel('gradient 1')
    % ylabel('gradient 2')
    % 
    % [c1,p1] = corr(dat,grad_cifti(1,:)','type','s') ; 
    % [c2,p2] = corr(dat,grad_cifti(2,:)','type','s') ; 
    % 
    % text(0.6,0.05,{ 'grad 1' [ 'rho: ' num2str(round(c1,2)) ', (' num2str(round(p1,2)) ')' ]},...
    %     'units','normalized') ; 
    % text(0.2,0.95,{ 'grad 2' [ 'rho: ' num2str(round(c2,2)) ', (' num2str(round(p2,2)) ')' ]},...
    %     'units','normalized') ; 
end

% set(gcf,'Position',[100 100 1000 600])
% set(gcf,'Color','w')
% orient(gcf,'landscape')

%%

% out_figdir = [ './reports/figures/figH/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/hub_pic.png' ] ; 
% print(filename,'-dpng')
% close(gcf)