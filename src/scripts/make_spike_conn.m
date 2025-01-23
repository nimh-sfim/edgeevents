%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_mats.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        spike_mats.(sdx{1}){idx} = readdat.spike_len_mat  ; 
    
    end

end

%% lets get some average mats for each spike len division

maxspk = 100 ; 

high_bin = 4 ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

spklen_names = {'short' 'inter' 'long'} ; 
spike_conn = struct() ; 

for sdx = subsets

    for jdx = spklen_names
        spike_conn.(sdx{1}).(jdx{1}) = zeros(finfo.nnodes) ; 
        spike_conn_indiv.(sdx{1}).(jdx{1}) = zeros(finfo.nnodes+55,finfo.nnodes+55,length(sublist.(sdx{1}))) ; 
    end

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        spkmat = spike_mats.(sdx{1}){idx} ; 

        dd = discretize(spike_mats.(sdx{1}){idx} ,lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ; 
        
        for ndx = 1:length(spklen_names) 
            nn = spklen_names{ndx} ; 

            tmpmat = mksq(count_spks(dd==ndx)) ; 

            spike_conn.(sdx{1}).(nn) = tmpmat(1:finfo.nnodes,1:finfo.nnodes) + ...
                spike_conn.(sdx{1}).(nn) ; 

            spike_conn_indiv.(sdx{1}).(nn)(:,:,idx) = tmpmat ; 

        end

    end   

    % divide by n subjects
    for ndx = 1:length(spklen_names) 
        nn = spklen_names{ndx} ; 

        spike_conn.(sdx{1}).(nn) = spike_conn.(sdx{1}).(nn) ./ length(sublist.(sdx{1})) ; 

    end

end

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

%% save the spike_conn

if strcmp(imglob,'REST1_LR')
    OUTSTR = 'sch200_2' ; 
end

filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
save(filename,'spike_conn','spklen_names','-v7.3')

filename = [ DD.PROC '/spk_conn_indiv_' OUTSTR '.mat' ] ; 
save(filename,'spike_conn_indiv','spklen_names','-v7.3')

filename = [ DD.PROC '/fc_conn_indiv_' OUTSTR '.mat' ] ; 
save(filename,'fc_conn_indiv','-v7.3')


%% the longer spike conn
% TODO: record each subject's super long individual matrix

spklen_names_ex = {'longer' 'longest'} ; 
spike_conn_ex =  struct() ; 
lll = [8 14] ;

spike_conn_ex_indiv = struct() ; 

for sdx = subsets

    for jdx = spklen_names_ex
        spike_conn_ex.(sdx{1}).(jdx{1}) = zeros(finfo.nnodes) ;
        spike_conn_ex_indiv.(sdx{1}).(jdx{1}) = zeros(255*254/2,length(sublist.(sdx{1}))) ; 
    end

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        spkmat = spike_mats.(sdx{1}){idx} ; 
        
        for ndx = 1:length(spklen_names_ex) 
            nn = spklen_names_ex{ndx} ; 

            tmpmat = mksq(count_spks(spkmat>=lll(ndx))) ; 

            spike_conn_ex.(sdx{1}).(nn) = tmpmat(1:finfo.nnodes,1:finfo.nnodes) + ...
                spike_conn_ex.(sdx{1}).(nn) ; 

            spike_conn_ex_indiv.(sdx{1}).(nn)(:,idx) = tv(tmpmat) ; 

        end

    end   

    % divide by n subjects
    for ndx = 1:length(spklen_names_ex) 
        nn = spklen_names_ex{ndx} ; 

        spike_conn_ex.(sdx{1}).(nn) = spike_conn_ex.(sdx{1}).(nn) ./ length(sublist.(sdx{1})) ; 

    end

end

filename = [ DD.PROC '/spk_conn_ex_avg_' OUTSTR '.mat' ] ; 
save(filename,'spike_conn_ex','spklen_names_ex','-v7.3')

filename = [ DD.PROC '/spk_conn_ex_indiv_' OUTSTR '.mat' ] ; 
save(filename,'spike_conn_ex_indiv','spklen_names_ex','-v7.3')

% run up to here for the sch200_2

%%

TL = tiledlayout(2,3) ;
TL.TileSpacing = 'compact'; 

spkken_names_long =  {'short' 'intermediate' 'long'} ; 

CM = flipud(gnbu(100)) ; 

for idx = 1:3

    nexttile()

    if idx == 1
        imsc_grid_comm(spike_conn.subset1.(spklen_names{idx}),parc.ca(1:finfo.nnodes), ...
            1,[1 1 1],[],parc.names(1:17))
        set(gca,'TickLength',[ 0 0])
    else
        imsc_grid_comm(spike_conn.subset1.(spklen_names{idx}),parc.ca(1:finfo.nnodes), ...
            1,[1 1 1],[],[])
        set(gca,'TickLength',[ 0 0])
        yticks('')
    end
    
    clim([0 15])

    axis square
    xticks('')
    xticklabels('')

    title(spkken_names_long{idx})

    if idx == 3
        cb = colorbar() ; 
        cb.Label.String = 'count' ; 
    end

end

% for idx = 1:3
% 
%     nexttile()
% 
%     if idx == 1
%         imsc_grid_comm(get_blocky(spike_conn.subset1.(spklen_names{idx}),...
%             parc.ca(1:finfo.nnodes)),1:17,1,[1 1 1],1, ...
%             arrayfun(@(i_) [ parc.names{i_} ': ' num2str(i_) ],1:17,'UniformOutput',false))
%         set(gca,'TickLength',[ 0 0])
%     else
%         imsc_grid_comm(get_blocky(spike_conn.subset1.(spklen_names{idx}),parc.ca(1:finfo.nnodes)),1:17,1,[1 1 1],1)
%         set(gca,'TickLength',[ 0 0])
% 
%         yticks(1:17)
%         yticklabels(cellstr(num2str((1:17)')))
% 
%     end
% 
%     clim([0 6])
% 
%     axis square
%     xticks(1:17)
%     xticklabels(cellstr(num2str((1:17)')))
% 
%     colormap(CM)
% 
% 
%     if idx == 3
%         cb = colorbar() ; 
%         cb.Label.String = 'count' ; 
%     end
% 
% end

%

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
remap_colors = remaplabs(1:17,g1sort,1:17) ; 
% cmap = get_nice_yeo_cmap('grad1') ; 
cmap = viridis(17) ; 

[u,v] = find(triu(ones(finfo.nnodes),1)) ; 
ca_u = remap_labs(u) ; ca_v = remap_labs(v) ; 

for idx = 1:3

    nexttile(TL)
    mat = (spike_conn.subset1.(spklen_names{idx})) ; 
    dat = [] ; 
    cat = [] ; 

    for jdx = 1:max(parc.ca(1:finfo.nnodes))
        mm = triu(mksq(logical((ca_u == jdx) & (ca_v == jdx))),1) ; 
        dat = [ dat ; mat(~~mm) ] ;
        cat = [ cat ; ones(length(mat(~~mm)),1).*jdx] ; 
    end

    % dat = mean(spike_conn.subset1.(spklen_names{idx}))' ; 

    % fcn_boxpts(dat(datplot_sinds),...
    %     datplot_ca,cmap(datplot_lab_sinds,:),...
    %     0,parc.names(g1sort))

    % fcn_boxpts(dat,...
    %     cat,cmap(datplot_lab_sinds,:),...
    %     0,parc.names(g1sort))

    fcn_boxpts(dat,...
        cat,cmap,...
        0,parc.names(g1sort))
    ylim([0 18])

    %ylim([ 0 .1])
    axis square

    set(gca, 'TickLabelInterpreter', 'none')


end

%%

set(gcf,'Position',[100 100 1000 600])

%%

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mats_separated.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% gather some stats

% for idx = 1:3
% 
%     mat = (spike_conn.subset1.(spklen_names{idx})) ; 
% 
%     corr(tv(mat),tv(meanfc),'type','s')
%     prctile(bootstrp(1000,@(x,y) corr(x,y,'type','s'),tv(mat),tv(meanfc)),[2.5 97.5])
% end

corr2fc.subset1 = nan(length(sublist.(sdx{1})),3) ; 
corr2fc.subset2 = nan(length(sublist.(sdx{1})),3) ; 

for sdx = subsets

    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)

       for ndx = 1:length(spklen_names) 
            nn = spklen_names{ndx} ;
    
            corr2fc.(sdx{1})(idx,ndx) = corr(tv(spike_conn_indiv.(sdx{1}).(nn)(:,:,idx)),...
                tv(fc_conn_indiv.(sdx{1})(:,:,idx)),'type','s')  ;
       end
    
    end

end

%%

tiledlayout(1,2)
cm = plasma(3) ; 

for sdx = subsets

    nexttile()

    for idx = 1:3
        histogram(corr2fc.(sdx{1})(:,idx),FaceColor=cm(idx,:))
        hold on
    end
    hold off

    axis square
    xlabel('rank corr. to FC')
    ylim([0 50])

    title(sdx{1})

    if strcmp(sdx{1},'subset2')
        legend(spklen_names,'Location','northeastoutside')
    else
        ylabel('num. subjects')
    end
end

set(gcf,'Position',[100 100 600 500])

%%

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mats_corr2fc_hist.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

filename = [ DD.PROC '/spk_mats_corr2fc_' OUTSTR '.mat' ] ; 
save(filename,'corr2fc','spklen_names','-v7.3')

%% and some permulation based anova on those mats


TL = tiledlayout(1,3) ;
TL.TileSpacing = 'compact'; 

spkken_names_long =  {'short' 'intermediate' 'long'} ; 

for idx = 1:3

    disp(idx)

    mat = (spike_conn.subset1.(spklen_names{idx})) ; 

    [aa,bb,cc] = perm_mat_anova1(mat,remap_labs,1000) ; 
    
    mm = multcompare(cc,'Alpha',0.001,'Display','off') ; 
    %mmm = mksq(mm(:,6)) ; 
    trilm = logical(tril(ones(17),-1)) ;
    mmm = zeros(17) ; 
    mmm(trilm) = mm(:,6) ;
    mmm = mmm + mmm' ; 

    nt = nexttile(TL)

    if idx == 1
        h = imsc_grid_comm(mmm,1:17,0.5,[1 1 1],[1 1 1],parc.names(g1sort)) ; 
    else
        h = imsc_grid_comm(mmm,1:17,0.5,[1 1 1],[1 1 1]) ; 
    end
    axis square
    colormap('parula')
    title(spkken_names_long{idx})

    sigthr = 1/(17*16/2) ; 
    h.AlphaData = (mmm <= sigthr)+.5 ; 

    xticks('')


    xlabel(['F(',num2str(round(aa.df_report(1),2,'sig')) ,', ',...
            num2str(aa.df_report(2)),') = ' num2str(aa.f) , ...
            ', p < ' , num2str(round(aa.p,2,'sig'))])

    if idx==3
        colorbar()
    end
end


%%

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[100 100 1000 600])

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_comms_anova_sig.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% lil gradient picutre

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

parc_plot_wcolorbar(grad_cifti(1,:),surfss,annotm,...
    [min(grad_cifti(1,:)) max(grad_cifti(1,:))],viridis(50),[100 100 600 1000])


set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/grad_1.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)


%%

CM = flipud(gnbu(100)) ; 

for idx = 1:3

    % [dat,~] = pca(spike_conn.subset1.(spklen_names{idx}),'NumComponents',1) ; 
    %[dat,~] = svds(spike_conn.subset1.(spklen_names{idx}),1) ;
    % dat = sum(spike_conn.subset1.(spklen_names{idx})) ; 
    % dat = eigenvector_centrality_und(spike_conn.subset1.(spklen_names{idx})) ; 
    dat = mean(spike_conn.subset1.(spklen_names{idx})) ; 

    parc_plot_wcolorbar(dat,surfss,annotm,...
        [min(dat) max(dat)],CM,[100 100 600 1000])

    set(gcf,'Color','w')

    
    out_figdir = [ './reports/figures/figA/' ]
    mkdir(out_figdir)
    filename = [out_figdir '/spike_cortex_' spklen_names{idx} '.png' ] ; 
    print(filename,'-dpng')
    close(gcf)

end

%% for supp

TL = tiledlayout(3,3)
TL.TileSpacing = 'compact'; 

spklen_names = {'short' 'inter' 'long'} ; 
spkken_names_long =  {'short' 'intermediate' 'long'} ; 

CM = flipud(gnbu(100)) ; 


for idx = 1:3

    nexttile()

    if idx == 1
        imsc_grid_comm(ranktrnsf_und(spike_conn.subset1.(spklen_names{idx})),parc.ca(1:finfo.nnodes), ...
            1,[1 1 1],[],parc.names(1:17))
        set(gca,'TickLength',[ 0 0])
    else
        imsc_grid_comm(ranktrnsf_und(spike_conn.subset1.(spklen_names{idx})),parc.ca(1:finfo.nnodes), ...
            1,[1 1 1],[],[])
        set(gca,'TickLength',[ 0 0])
        yticks('')
    end
   
    axis square
    xticks('')
    xticklabels('')

    title(spkken_names_long{idx})

    colormap(CM)

    if idx == 3
        cb = colorbar() ; 
        cb.Label.String = 'count' ; 
    end

end

% and compare the mats to each others
for idx = 1:3
    for jdx = idx+1:3
    
        nexttile()

        name1 = spkken_names_long{idx} ; 
        name2 = spkken_names_long{jdx} ; 

        dat1 = spike_conn.subset1.(spklen_names{idx}) ; 
        dat2 = spike_conn.subset1.(spklen_names{jdx}) ; 

        ss = scatter(tv(dat1),tv(dat2),"filled",'MarkerEdgeColor','none',...
            'MarkerFaceColor',[0.8 0.8 0.8]) ; 
        rho = corr(tv(dat1),tv(dat2),'Type','Spearman') ; 
        text(0.05,0.9,['rho: ' num2str(round(rho,2))],...
            'Units','normalized')

        xlabel(name1) 
        ylabel(name2)

        rl = refline(1,0)
        rl.LineWidth = 2 ;

    end
end

for idx = 1:3
    
    nexttile()

    name1 = spkken_names_long{idx} ; 

    dat1 = spike_conn.subset1.(spklen_names{idx}) ; 
    dat2 = spike_conn.subset2.(spklen_names{idx}) ; 

    ss = scatter(tv(dat1),tv(dat2),"filled",'MarkerEdgeColor','none',...
        'MarkerFaceColor',[0.8 0.8 0.8]) ; 
    rho = corr(tv(dat1),tv(dat2),'Type','Spearman') ; 
    text(0.05,0.9,['rho: ' num2str(round(rho,2))],...
        'Units','normalized')

    title(name1) 
    xlabel('subset 1')
    ylabel('subset 2')

    rl = refline(1,0)
    rl.LineWidth = 2 ;

end

set(gcf,'Position',[100 100 1200 1000])

set(gcf,'Color','w')
orient(gcf,'landscape')

%%

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mats_rank_n_scatter.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

% %% compare the short and long spike spaces
% 
% TL = tiledlayout(1,2)
% 
% nexttile()
% 
% long_degree = sum(spike_conn.subset1.long) ; 
% short_degree = sum(spike_conn.subset1.short) ; 
% 
% tr_degree = [tiedrank(long_degree(:))  tiedrank(short_degree(:)) ] ; 
% [~,ii] = max(tr_degree,[],2) ; 
% 
% dtl_degree = arrayfun(@(i_) point_to_line([tr_degree(i_,1) tr_degree(i_,2) 0],[0 0 0],[1 1 0]),1:length(tr_degree)) ; 
% dtl_degree(ii==2) = -dtl_degree(ii==2) ; 
% 
% 
% scatter(tr_degree(:,1),tr_degree(:,2),40,dtl_degree,'filled','MarkerEdgeColor','flat')
% colormap(rdbu)
% clim([-max(abs(dtl_degree)) max(abs(dtl_degree))])
% rl = refline(1,0)
% rl.Color = [0.8 0.8 0.8 ] ; 
% rl.LineWidth = 2 ;
% 
% axis square
% 
% xlabel('long degree rank')
% ylabel('short degree rank')
% 
% colorbar
% 
% nexttile()
% 
% long_edge = tiedrank(tv(spike_conn.subset1.long)) ; 
% short_edge = tiedrank(tv(spike_conn.subset1.short)) ; 
% 
% tr_edge = [tiedrank(long_edge(:))  tiedrank(short_edge(:)) ] ; 
% [~,ii] = max(tr_edge,[],2) ; 
% 
% dtl_edge = arrayfun(@(i_) point_to_line([tr_edge(i_,1) tr_edge(i_,2) 0],[0 0 0],[1 1 0]),1:length(tr_edge)) ; 
% dtl_edge(ii==2) = -dtl_edge(ii==2) ; 
% 
% dtl_edge_mat = mksq(dtl_edge) ; 
% 
% scatter(tr_edge(:,1),tr_edge(:,2),6,dtl_edge,'filled','MarkerEdgeColor','flat')
% rl = refline(1,0)
% rl.Color = [0.8 0.8 0.8 ] ; 
% rl.LineWidth = 2 ;
% 
% clim([-max(abs(dtl_edge)) max(abs(dtl_edge))])
% 
% axis square 
% 
% xlabel('long edge rank')
% ylabel('short edge rank')
% 
% colorbar
% 
% set(gcf,'Position',[100 100 800 500])
% 
% %%
% 
% set(gcf,'Color','w')
% orient(gcf,'landscape')
% 
% out_figdir = [ './reports/figures/figA/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/shortlong_pref_scatter.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)
% 
% %% and make the cortex plot
% 
% TL = tiledlayout(2,1)
% nt1 = nexttile()
% 
% pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dtl_degree ,...
%     'valRange',[-max(abs(dtl_degree)) max(abs(dtl_degree))],...
%     'cmap',rdbu, ...
%     'viewcMap',0,'newFig',0,'viewStr','all',...
%     'parenth',TL)
% pp.Layout = nt1.Layout ; 
% 
% nt2 = nexttile(TL)
% 
% hh = imagesc(dtl_degree) 
% cb = colorbar()
% clim([-max(abs(dtl_degree)) max(abs(dtl_degree))])
% hh.Visible = 'off' ;
% hh.Parent.Visible = 'off' ; 
% cb.Location = "north" ; 
% cm = colormap() ; 
% colormap(nt2,cm(2:end,:))
% 
% TL.TileSpacing = 'tight'
% 
% set(gcf,'Position',[100 100 600 1000])
% 
% out_figdir = [ './reports/figures/figA/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/shortlong_pref_cortex.png' ] ; 
% print(filename,'-dpng')
% close(gcf)

% %% top long edges (depreciated)
% 
% pp = 5 ; 
% 
% topedges_mat = dtl_edge_mat>prctile(dtl_edge_mat(:),100-pp) ; 
% 
% TL = tiledlayout(1,2)
% 
% nexttile()
% 
% h = imsc_grid_comm(topedges_mat,parc.ca(1:finfo.nnodes), ...
%     1,[1 1 1],[],parc.names(1:17))
% set(gca,'TickLength',[ 0 0])
% 
% clim([0 1])
% axis square
% 
% b = blues(11) ; 
% bb = b(6,:) ; 
% 
% colormap(gca,[1 1 1 ; bb])
% 
% xticks('')
% 
% nexttile()
% 
% imsc_grid_comm(get_blocky(topedges_mat,parc.ca(1:finfo.nnodes)),1:17, ...
%     1,[1 1 1],[],parc.names(1:17))
% set(gca,'TickLength',[ 0 0])
% 
% colormap(gca,interp_cmap([1 1 1], bb,50))
% 
% axis square
% 
% xticks(1:17)
% 
% TL.Title.String = { [ 'Top ' num2str(pp) '% long-perference edges' ]  ''  } 
% 
% set(gcf,'Position',[100 100 800 500])
% set(gcf,'Color','w')
% orient(gcf,'landscape')
% 
% %%
% 
% out_figdir = [ './reports/figures/figA/' ] 
% mkdir(out_figdir)
% filename = [out_figdir '/top_longpref_edges.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)
% 
% 
% %%
% % 
% % figure
% % imsc_grid_comm(get_blocky(dtl_edge_mat,parc.ca(1:finfo.nnodes)),1:17, ...
% %     1,[1 1 1],[],parc.names(1:17))
% % set(gca,'TickLength',[ 0 0])
% % 
% % clim([-max(abs(dtl_edge)) max(abs(dtl_edge))])
% % 
% % axis square
% % 
% % %% 
% % 
% % figure
% % imsc_grid_comm(dtl_edge_mat>prctile(dtl_edge_mat(:),95),parc.ca(1:finfo.nnodes), ...
% %     1,[1 1 1],[],parc.names(1:17))
% % set(gca,'TickLength',[ 0 0])
% % 
% % clim([0 1])
% % 
% % axis square
% % 
% % %%
% 
% %% compare to the two gradients! 
% 
% 
% grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 
% 
% scatter(grad_cifti(1,:),grad_cifti(2,:),40,dtl_degree,'filled','MarkerEdgeColor','flat')
% colormap(rdbu)
% 
% xlabel('gradient 1')
% ylabel('gradient 2')
% 
% [r1,p1] = corr(grad_cifti(1,:)',dtl_degree(:),'type','s') ; 
% [r2,p2] = corr(grad_cifti(2,:)',dtl_degree(:),'type','s') ; 
% 
% text(0.6,0.9,[ 'rho w/ grad. 1: ' num2str(round(r1,2)) ],'Units','normalized')
% text(0.6,0.85,[ 'rho w/ grad. 2: ' num2str(round(r2,2)) ],'Units','normalized')
% 
% axis square
% 
% set(gcf,'Position',[100 100 400 400])
% 
% %%
% 
% out_figdir = [ './reports/figures/figA/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/shortlong_pref_grad12.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)
% 
% %%
% 
% grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 
% 
% dat = spike_conn.subset1.long ./ (spike_conn.subset1.short + spike_conn.subset1.inter + spike_conn.subset1.long) ; 
% dat(isnan(dat)) = 0 ;
% dat = sum(dat) ; 
% 
% scatter3(grad_cifti(1,:),grad_cifti(2,:),grad_cifti(3,:),40, dat ,'filled','MarkerEdgeColor','flat')
% 
% xlabel('gradient 1')
% ylabel('gradient 2')
% 
% clim([ prctile(dat,[02.5])  prctile(dat,97.5) ])
% 
% % [r1,p1] = corr(grad_cifti(1,:)',dtl_degree(:),'type','s') ; 
% % [r2,p2] = corr(grad_cifti(2,:)',dtl_degree(:),'type','s') ; 
% % 
% % text(0.6,0.9,[ 'rho w/ grad. 1: ' num2str(round(r1,2)) ],'Units','normalized')
% % text(0.6,0.85,[ 'rho w/ grad. 2: ' num2str(round(r2,2)) ],'Units','normalized')
% 
% axis square
% 
% colormap(parula)
% 
% set(gcf,'Position',[100 100 400 400])
% 
% %% do a scatter to many more gradients
% 
% tiledlayout(3,3)
% 
% for idx = 1:9
% 
%     nexttile()
% 
%     scatter_w_rho(grad_cifti(idx,:)',dtl_degree(:),"filled",'MarkerEdgeColor','none',...
%             'MarkerFaceColor',[0.8 0.8 0.8])
%     axis square
% 
%     title({['gradient ' num2str(idx) ' vs'] 'long-short pref' })
% 
% 
%     if idx == 4
%         ylabel('long-short preference')
%     end
% 
%     if idx == 8
%         xlabel('gradient val (a.u.)')
%     end
% 
% end 
% 
% set(gcf,'Position',[100 100 600 600])
% 
% out_figdir = [ './reports/figures/figA/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/shortlong_pref_allgrad.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)

%% just some random stuff

% res = zeros(199,1200) ; 
% resfc = zeros(200,1) ; 
% 
% for idx = 1:100
% 
%     disp(idx)
% 
%     tmpfc = corr(datStr(idx).ts) ; 
%     resfc = resfc + tmpfc(1:200,142) ; 
% 
%     tmpets = get_ets(datStr(idx).ts) ; 
%     res = tmpets(:,get_ets_node_inds([1:200] == 142))' + res ; 
% 
% end
% 
% resfc = resfc./100 ; 
% resfc(142) = [] ; 
% 
% res = res ./ 100 ; 
% 
% dat = ets(:,get_ets_node_inds([1:200] == 142)) ;
% 
% [uu,vv] = pca(dat','NumComponents', 1) ; 
% 
% figure , imagesc(dat(sortedInd(uu),:)') ; colormap(rdbu) ; clim([-8 8]) ; colorbar

