%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2.25 ; 

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_mats.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
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

        end

    end   

    % divide by n subjects
    for ndx = 1:length(spklen_names) 
        nn = spklen_names{ndx} ; 

        spike_conn.(sdx{1}).(nn) = spike_conn.(sdx{1}).(nn) ./ length(sublist.(sdx{1})) ; 

    end

end

%% save the spike_conn

filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
save(filename,'spike_conn','spklen_names','-v7.3')

%%

TL = tiledlayout(2,3) ;
TL.TileSpacing = 'compact'; 

spkken_names_long =  {'short' 'intermediate' 'long'} ; 

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
    
    clim([0 8])

    axis square
    xticks('')
    xticklabels('')

    title(spkken_names_long{idx})

    colormap(flipud(inferno()))

    if idx == 3
        cb = colorbar() ; 
        cb.Label.String = 'count' ; 
    end

end

for idx = 1:3

    nexttile()

    if idx == 1
        imsc_grid_comm(get_blocky(spike_conn.subset1.(spklen_names{idx}),...
            parc.ca(1:finfo.nnodes)),1:17,1,[1 1 1],1, ...
            arrayfun(@(i_) [ parc.names{i_} ': ' num2str(i_) ],1:17,'UniformOutput',false))
        set(gca,'TickLength',[ 0 0])
    else
        imsc_grid_comm(get_blocky(spike_conn.subset1.(spklen_names{idx}),parc.ca(1:finfo.nnodes)),1:17,1,[1 1 1],1)
        set(gca,'TickLength',[ 0 0])

        yticks(1:17)
        yticklabels(cellstr(num2str((1:17)')))

    end
    
    clim([0 8])

    axis square
    xticks(1:17)
    xticklabels(cellstr(num2str((1:17)')))
    
    % colormap(flipud(inferno()))
    colormap(flipud(viridis()))


    if idx == 3
        cb = colorbar() ; 
        cb.Label.String = 'count' ; 
    end

end

set(gcf,'Position',[100 100 800 500])

%%

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mats_separated.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% for supp

TL = tiledlayout(2,3)
TL.TileSpacing = 'compact'; 

spklen_names = {'short' 'inter' 'long'} ; 
spkken_names_long =  {'short' 'intermediate' 'long'} ; 

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

    colormap(flipud(inferno()))

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

%%

set(gcf,'Position',[100 100 800 500])

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_mats_rank_n_scatter.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% compare the short and long spike spaces

TL = tiledlayout(1,2)

nexttile()

long_degree = sum(spike_conn.subset1.long) ; 
short_degree = sum(spike_conn.subset1.short) ; 

tr_degree = [tiedrank(long_degree(:))  tiedrank(short_degree(:)) ] ; 
[~,ii] = max(tr_degree,[],2) ; 

dtl_degree = arrayfun(@(i_) point_to_line([tr_degree(i_,1) tr_degree(i_,2) 0],[0 0 0],[1 1 0]),1:length(tr_degree)) ; 
dtl_degree(ii==2) = -dtl_degree(ii==2) ; 


scatter(tr_degree(:,1),tr_degree(:,2),40,dtl_degree,'filled','MarkerEdgeColor','flat')
colormap(rdbu)
clim([-max(abs(dtl_degree)) max(abs(dtl_degree))])
rl = refline(1,0)
rl.Color = [0.8 0.8 0.8 ] ; 
rl.LineWidth = 2 ;

axis square

xlabel('long degree rank')
ylabel('short degree rank')

colorbar

nexttile()

long_edge = tiedrank(tv(spike_conn.subset1.long)) ; 
short_edge = tiedrank(tv(spike_conn.subset1.short)) ; 

tr_edge = [tiedrank(long_edge(:))  tiedrank(short_edge(:)) ] ; 
[~,ii] = max(tr_edge,[],2) ; 

dtl_edge = arrayfun(@(i_) point_to_line([tr_edge(i_,1) tr_edge(i_,2) 0],[0 0 0],[1 1 0]),1:length(tr_edge)) ; 
dtl_edge(ii==2) = -dtl_edge(ii==2) ; 

dtl_edge_mat = mksq(dtl_edge) ; 

scatter(tr_edge(:,1),tr_edge(:,2),6,dtl_edge,'filled','MarkerEdgeColor','flat')
rl = refline(1,0)
rl.Color = [0.8 0.8 0.8 ] ; 
rl.LineWidth = 2 ;

clim([-max(abs(dtl_edge)) max(abs(dtl_edge))])

axis square 

xlabel('long edge rank')
ylabel('short edge rank')

colorbar

set(gcf,'Position',[100 100 800 500])

%%

set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/shortlong_pref_scatter.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% and make the cortex plot

TL = tiledlayout(2,1)
nt1 = nexttile()

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dtl_degree ,...
    'valRange',[-max(abs(dtl_degree)) max(abs(dtl_degree))],...
    'cmap',rdbu, ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL)
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL)

hh = imagesc(dtl_degree) 
cb = colorbar()
clim([-max(abs(dtl_degree)) max(abs(dtl_degree))])
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight'

set(gcf,'Position',[100 100 600 1000])

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/shortlong_pref_cortex.png' ] ; 
print(filename,'-dpng')
close(gcf)

%% top long edges 

pp = 5 ; 

topedges_mat = dtl_edge_mat>prctile(dtl_edge_mat(:),100-pp) ; 

TL = tiledlayout(1,2)

nexttile()

h = imsc_grid_comm(topedges_mat,parc.ca(1:finfo.nnodes), ...
    1,[1 1 1],[],parc.names(1:17))
set(gca,'TickLength',[ 0 0])

clim([0 1])
axis square

b = blues(11) ; 
bb = b(6,:) ; 

colormap(gca,[1 1 1 ; bb])

xticks('')

nexttile()

imsc_grid_comm(get_blocky(topedges_mat,parc.ca(1:finfo.nnodes)),1:17, ...
    1,[1 1 1],[],parc.names(1:17))
set(gca,'TickLength',[ 0 0])

colormap(gca,interp_cmap([1 1 1], bb,50))

axis square

xticks(1:17)

TL.Title.String = { [ 'Top ' num2str(pp) '% long-perference edges' ]  ''  } 

set(gcf,'Position',[100 100 800 500])
set(gcf,'Color','w')
orient(gcf,'landscape')

%%

out_figdir = [ './reports/figures/figA/' ] 
mkdir(out_figdir)
filename = [out_figdir '/top_longpref_edges.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)


%%
% 
% figure
% imsc_grid_comm(get_blocky(dtl_edge_mat,parc.ca(1:finfo.nnodes)),1:17, ...
%     1,[1 1 1],[],parc.names(1:17))
% set(gca,'TickLength',[ 0 0])
% 
% clim([-max(abs(dtl_edge)) max(abs(dtl_edge))])
% 
% axis square
% 
% %% 
% 
% figure
% imsc_grid_comm(dtl_edge_mat>prctile(dtl_edge_mat(:),95),parc.ca(1:finfo.nnodes), ...
%     1,[1 1 1],[],parc.names(1:17))
% set(gca,'TickLength',[ 0 0])
% 
% clim([0 1])
% 
% axis square
% 
% %%

%% compare to the two gradients! 


grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

scatter(grad_cifti(1,:),grad_cifti(2,:),40,dtl_degree,'filled','MarkerEdgeColor','flat')
colormap(rdbu)

xlabel('gradient 1')
ylabel('gradient 2')

[r1,p1] = corr(grad_cifti(1,:)',dtl_degree(:),'type','s') ; 
[r2,p2] = corr(grad_cifti(2,:)',dtl_degree(:),'type','s') ; 

text(0.6,0.9,[ 'rho w/ grad. 1: ' num2str(round(r1,2)) ],'Units','normalized')
text(0.6,0.85,[ 'rho w/ grad. 2: ' num2str(round(r2,2)) ],'Units','normalized')

axis square

set(gcf,'Position',[100 100 400 400])

%%

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/shortlong_pref_grad12.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%

grad_cifti = squeeze(niftiread('data/external/hpc_grad_sch200-yeo17.pscalar.nii')) ; 

dat = spike_conn.subset1.long ./ (spike_conn.subset1.short + spike_conn.subset1.inter + spike_conn.subset1.long) ; 
dat(isnan(dat)) = 0 ;
dat = sum(dat) ; 

scatter3(grad_cifti(1,:),grad_cifti(2,:),grad_cifti(3,:),40, dat ,'filled','MarkerEdgeColor','flat')

xlabel('gradient 1')
ylabel('gradient 2')

clim([ prctile(dat,[02.5])  prctile(dat,97.5) ])

% [r1,p1] = corr(grad_cifti(1,:)',dtl_degree(:),'type','s') ; 
% [r2,p2] = corr(grad_cifti(2,:)',dtl_degree(:),'type','s') ; 
% 
% text(0.6,0.9,[ 'rho w/ grad. 1: ' num2str(round(r1,2)) ],'Units','normalized')
% text(0.6,0.85,[ 'rho w/ grad. 2: ' num2str(round(r2,2)) ],'Units','normalized')

axis square

colormap(parula)

set(gcf,'Position',[100 100 400 400])

%% do a scatter to many more gradients

tiledlayout(3,3)

for idx = 1:9

    nexttile()

    scatter_w_rho(grad_cifti(idx,:)',dtl_degree(:),"filled",'MarkerEdgeColor','none',...
            'MarkerFaceColor',[0.8 0.8 0.8])
    axis square

    title({['gradient ' num2str(idx) ' vs'] 'long-short pref' })


    if idx == 4
        ylabel('long-short preference')
    end

    if idx == 8
        xlabel('gradient val (a.u.)')
    end

end 

set(gcf,'Position',[100 100 600 600])

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/shortlong_pref_allgrad.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

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
