%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

addpath('/Users/faskowitzji/Documents/MATLAB/conn/')
addpath('/Users/faskowitzji/Documents/MATLAB/spm12/')

SPK_THR = 2.25 ; 

%%

for idx = 1:NSUBS

    disp(idx)

    filename = [DD.PROC '/' imglob '/' datStr(idx).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 

    if isfile(filename)
        disp(['already finsiehd sub: ' datStr(idx).sub ])
        continue
    end

    tmpts = datStr(idx).ts ; 
    tmpets = get_ets(tmpts) ; 

    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ; 

    mkdir([DD.PROC '/' imglob '/'])
    save(filename,...
        'spike_len_cell','spike_len_mat','-v7.3')

end

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_lengths.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
    
        spike_lengths.(sdx{1}){idx} = readdat.spike_len_cell  ; 
    
    end

end

for sdx = subsets
    tmp = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.(sdx{1}){i_}')),1:length(spike_lengths.(sdx{1}))','UniformOutput',false)') ;
    % get rid of 0 lengths
    spike_lengths.all.(sdx{1}) = nonzeros(tmp) ; 
end

%% look at spike length for each pers

% spike_prct.subset2.prct = nan(100,length(sublist.subset2)) ; 
spike_prct = struct() ;  
spike_prct.all.subset1 = prctile(spike_lengths.all.subset1, 1:100 ) ;
spike_prct.all.subset2 = prctile(spike_lengths.all.subset2, 1:100 ) ; 

for sdx = subsets

    spike_prct.(sdx{1}).prct = nan(100,length(sublist.(sdx{1}))) ; 

    for idx = 1:length(sublist.(sdx{1}))
        disp(idx)
        tmp = nonzeros(single(cell2mat(spike_lengths.(sdx{1}){idx}'))) ; 
        spike_prct.(sdx{1}).prct(:,idx) = prctile(tmp, 1:100 )  ; 
    end

end

% %% 
% 
% % look at the subjects 
% imagesc(spike_prct.subset1.prct(:,fliplr(sortedInd(mean(spike_prct.subset1.prct(95:100,:))))))

%% find way to naturally bin it
% look at this for fun, not a super intuitive figure

plot(mean(spike_prct.subset1.prct,2),'LineWidth',2)
rr = round(mean(spike_prct.subset1.prct,2)) ; 
urr = unique(rr) ;
nums = 1:100 ; 
for idx = 1:length(unique(rr))
    hold on
    ll = rr==urr(idx) ; 
    plot(nums(ll),ones(sum(ll),1).*urr(idx),...
        Color=[0.5 0.5 0.5 0.5],LineWidth=5)
end
hold off
ylim([0 15])

%% bin it in a data driven way

maxspk = 100 ; 

spike_dist = discretize(spike_lengths.all.subset1,1:1:maxspk) ; 
spike_tab = tabulate(spike_dist) ; 

% how much in lowest spike
low_bin_pct = spike_tab(1,3) ; 

% now find a high bin to get as near as possible
[~,mi] = min(abs(flipud(cumsum(flipud(spike_tab(2:end,3))))-low_bin_pct)) ; 
high_bin = mi + 1 ; % add 1 because we were aready looking at 2:end

% the low-med-high bins
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

%% 

% tiledlayout(1,2)

cm = plasma(3) ; 

clf

dd = discretize(spike_lengths.all.subset1,lowmedhigh_edges) ; 

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:40,FaceColor=cm(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    sum(dd==idx)./length(dd)
    hold on

end
hold off

xlim([1 20.5])

xlabel('spike length (sec)')
ylabel('count')

xticks((1:20)+0.5)
xticklabels(cellstr(num2str((1:20)'.*finfo.TR)))
xtickangle(45)

axis square

hist_gca = gca ;
group_lab = { ...
    strcat('short (',num2str(round(sum(dd==1)/length(dd)*100)),'%)') ... 
    strcat('intermed. (',num2str(round(sum(dd==2)/length(dd)*100)),'%)') ... 
    strcat('long (',num2str(round(sum(dd==3)/length(dd)*100)),'%)') ... 
    }

legend(hist_gca,group_lab,'Location','southeast')

aa = axes('Parent',gcf,'Position',[.48 .5 .4 .38],'Units','normalized')
box on
plot(cumsum(spike_tab(1:20,3)),'.-','MarkerSize',15,'Color',[0.5 0.5 0.5])
xlabel('spike length (sec)')
ylabel('cumulative percent')
xlim([1 20])
ylim([1 100])
axis square
xticks([1 5 10 15 20])
xticklabels(cellstr(num2str(([1 5 10 15 20])'.*finfo.TR)))

yticks([0 20 40 60 80 100])


set(gcf,'Position',[100 100 400 400])

%% 

set(gcf,'Color','w')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_hist.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% lets also make a log version of that plot

clf

dd = discretize(spike_lengths.all.subset1,lowmedhigh_edges) ; 

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:40,FaceColor=cm(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    hold on

end
hold off
set(gca,'YScale','log')

endval=40 ;
xlim([1 endval+0.5])

xlabel('spike length (sec)')
ylabel('count')

xticks((1:endval)+0.5)
xticklabels(cellstr(num2str((1:endval)'.*finfo.TR)))
xtickangle(45)

set(gcf,'Position',[100 100 800 800])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_hist_log.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)

%% lets do reproducilbility

clf

TL1 = tiledlayout(1,3)

for tdx = 1:2

    nexttile(TL1)
    
    TL2 = tiledlayout(TL1,2,1,'TileSpacing','none') ;

    TL2.Layout.Tile = tdx ; 
    TL2.Layout.TileSpan = [ 1 1] ;

    axis off 
    ss = {'subset1' 'subset2'} ; 
    
    for sdx = 1:2
        nt = nexttile(TL2) ; 
        dd = discretize(spike_lengths.all.(ss{sdx}),lowmedhigh_edges) ; 
        for idx = 1:max(dd)
        
            histogram(spike_lengths.all.(ss{sdx})(dd==idx), 0:1:40,FaceColor=cm(idx,:)) ;
            % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
            hold on
        
        end
        hold off
    
        if sdx ==2
        set(gca, 'ydir', 'reverse')
            xticks(1:2:40) 
            xtickangle(45)
            xx = xticklabels() ; 
            xticklabels(cellstr(num2str((str2double(xx).*finfo.TR))))
        else
            xticklabels([])
            title('histogram from each subset')
        end
    
        ylabel([ ss{sdx} ])
    
        if tdx == 2
            set(gca,'YScale','log')

        end

    end

end

% TL2.Title.String = 'histogram counts btwn. subsets'


nexttile(TL1)


p1 = discretize(spike_lengths.all.subset1,1:40) ; 
p2 = discretize(spike_lengths.all.subset2,1:40) ; 
% t1 = tabulate(p1) ; 
% t2 = tabulate(p2) ; 
ss = zeros(40,2) ; 
for idx = 1:40
   ss(idx,1) = sum(p1==idx) ; 
   ss(idx,2) = sum(p2==idx) ;
end

s = scatter(ss(:,1),ss(:,2),30,1:40,'filled') ; 
set(gca,'YScale','log')
set(gca,'XScale','log')
cb = colorbar() ;
colormap(parula(40))
cb.Label.String = 'event length' ;

cb.Ticks = [1:5:40] ; 
ll = str2double(cb.TickLabels()) ; 
cb.TickLabels = (cellstr(num2str(ll.*finfo.TR))) ;

set(gcf,'Color','w')

xlabel('subset 1')
ylabel('subset 2')

refline(1,0)

title('bin counts btwn. subsets')

%%

set(gcf,'Position',[100 100 1000 500])
set(gcf,'Color','w')
orient('landscape')

out_figdir = [ './reports/figures/figA/' ]
mkdir(out_figdir)
filename = [out_figdir '/spike_hist_subsets.pdf' ] ; 
print(filename,'-dpdf','-vector','-bestfit')
close(gcf)


%%

testmat = zeros(200) ; 
testmat(1,3) = 1 ;

testmat(1,101) = 1 ;
testmat(2,50) = 1 ;


pp = parula(3) ; 

rois.sph_r = ones(200,1); 
viz_conn_glassbrain(testmat,pp,rois)

%% lets look at the null dataaaa

% filename = [DD.PROC '/surrogate_' OUTSTR '_' , num2str(SPK_THR) , '_spk.mat'] ; 
% load(filename)

%% power-law ish??



