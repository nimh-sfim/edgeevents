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

% fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

refparc = parc.ca(1:200) ; 
nedges = ( finfo.nnodes * (finfo.nnodes-1) ) / 2 ; 

TR = 0.72 ; 

%% multiply two hrfs

% scratch
% hrf1 = getcanonicalhrf(3,1) ; 
% hrf2 = getcanonicalhrf(5,1) ; 
% 
% 
% rng(42)
% [~,randts] = pad_conv(rand(10000,2) > 0.7,hrf1) ; 
% 
% [~,spkc1] = spk_lenmat(randts(:,1) > prctile(randts(:,1),80)) ; 
% [~,spkc2] = spk_lenmat(randts(:,2) > prctile(randts(:,2),80)) ; 
% 
% 
% randets = get_ets(randts) ; 
% [~,etsspkc] = spk_lenmat(randets > prctile(randets,80)) ; 

%%

tiledlayout(2,1)

spklens = 1:2:29 ; 
CM = internet(max(spklens)) ; 

nexttile

for ii = spklens

    hh = getcanonicalhrf(ii,TR) ; 
    % ll = length(hh) ; 
    % % do a fair zscore, relfecting baseline at '0'
    % tmp = zscore([ hh -hh]) ; 
    % hh = tmp(1:ll) ; 

    plot((0:length(hh)-1).*TR,hh,'Color',CM(ii,:),'LineWidth',2)
    title('hrf')
    xlim([0 45])
    hold on

end
hold off

cb = colorbar()
colormap(CM) 
numcolors = length(spklens) ; 
cb.Ticks = (0+0.5*(numcolors)/numcolors:(numcolors)/numcolors:numcolors) ./ numcolors ;
cb.TickLabels = cellstr(num2str((spklens)'))
cb.Label.String = 'simulus duration (sec)'

ylabel('normalized amplitude')

%cb = colorbar() ; 

nexttile
for ii = spklens

    hh = getcanonicalhrf(ii,TR) ; 
    ll = length(hh) ; 
    % do a fair zscore, relfecting baseline at '0'
    tmp = zscore([ hh -hh]) ; 
    hh = tmp(1:ll) ; 

    hhh = prod([hh' hh'],2) ;

    plot((0:length(hhh)-1).*0.72,hhh,'Color',CM(ii,:),'LineWidth',2)
    xlim([0 45])
    hold on

end
hold off
title('hrf x hrf')

xlabel('time (sec)')
ylabel('co-fluctuation')

set(gcf,'Position',[100 100 800 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrf_hrfxhrf.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

tiledlayout(2,1)


hrf1 = getcanonicalhrf(1,TR) ; 
ll = length(hrf1) ; 
% do a fair zscore, relfecting baseline at '0'
tmp = zscore([ hrf1 -hrf1]) ; 
hrf1 = tmp(1:ll) ; 

nexttile()

plot((0:ll-1).*TR,getcanonicalhrf(1,TR),'LineWidth',2)
xlim([0 30])

title('hrf 1 sec. stim.')

nexttile()

plot((0:ll-1).*TR,prod([ hrf1' hrf1' ],2),'LineWidth',2)

for idx = 0:10

    %hh1 = [  hrf1'  hrf1' ] ; 
    hh2 = [ [ hrf1' ; zeros(idx,1)] [ zeros(idx,1) ; hrf1'] ] ; 
    tmp = zscore( [hh2 ; -hh2 ]) ; 
    hh2 = tmp(1:ll,:) ; 
    

    hold on
    plot((0:ll-1).*TR,prod(hh2,2),'LineWidth',2,'Color',CM(idx+5,:))
    
    title('hrf x hrf')
    xlim([0 30])

end
hold off

cb = colorbar()
colormap(CM(5:15,:)) 
numcolors = 11 ; 
cb.Ticks = (0+0.5*(numcolors)/numcolors:(numcolors)/numcolors:numcolors) ./ numcolors ;
cb.TickLabels = cellstr(num2str(((0:10).*TR)'))
cb.Label.String = 'offsec sec'


xlabel('time (sec)')

set(gcf,'Position',[100 100 800 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrf_hrfxhrf_offset.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% lets look at hrf values with differing starting points

TL = tiledlayout(3,3,'TileIndexing','columnmajor') ; 

rng(42)

nn = 25 ; 

%randoffsets = (rand(nn,2)-0.5) .* 5 ;
randoffsets = normrnd(0,1,[nn 2]) ; 

hrf1 = getcanonicalhrf(1,0.72) ; 

for shiftn = 1:3

    nexttile([1 2])
    
    for idx = 1:nn
    
        t1 = [zeros(shiftn,1) ; hrf1' ] +randoffsets(idx,1) ; 
        t2 = [ hrf1' ; zeros(shiftn,1) ] +randoffsets(idx,2) ; 
        ee = prod([t1 t2],2) ; 

        plot((0:length(ee)-1).*TR,ee,'Color',[ CM(20,:) 0.5],'LineWidth',2)
        hold on
        
    end
    hold off
    xlim([0 20])
    yline(2.25,'Color','b','LineWidth',1)
    
    title(['hrf x hrf with random start mag., offset ' num2str(shiftn.*TR) ' sec.'])

    if shiftn == 2
        xlabel('time (sec)')
    end

end

% lets do a simulation for each
nsim = 5000 ; 
spklen = zeros(nsim,1) ; 
%randoffsets = (rand(nsim,2)-0.5) .* 4 ; 
randoffsets = normrnd(0,1,[nsim 2]) ; 

histedges = [0:1:50].*TR ; 

for shiftn = 1:3

    nexttile
    
    for idx = 1:nsim
    
        t1 = [zeros(shiftn,1) ; hrf1' ] +randoffsets(idx,1) ; 
        t2 = [ hrf1' ; zeros(shiftn,1) ] +randoffsets(idx,2) ; 
        ee = prod([t1 t2],2) ; 
        %ee = ee(10:30) ; % only focus on hrf rise
    
        [~,ss ] = spk_lenmat(ee>2.25) ; 
        ss=cell2mat(ss) ; 

        % if length(ss) > 1
        %     disp('feff')
        %     ss
        % end

        if ~isnan(ss)
            spklen(idx) =  max(ss) ; % max, in case theres more than 1 for some reason 
        end
    end
    
    % there are bound to be some weird trials where above 0 the whole
    % time... just keep them in, they don't affect much
    % spklen(snpkle>20) = nan ;

    spklen = spklen .* TR ;
    nt = sum(spklen>0) ;

    histogram(nonzeros(spklen(spklen<(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]) ; 
    hold on
    histogram(nonzeros(spklen(spklen>=(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]) ; 
    hold off

    xticks([5 10 15])
    xticklabels({ '5' '10' '15' })

    xlim([1 20])

    % comput hitrate >= 2.88 sec
    text(0.8,0.8,{ 'events >= 2.88' [ num2str(round(sum(spklen>=(4.*TR)) ./ nt,2)) '%']} , ...
        'units','normalized','HorizontalAlignment','center') ; 

    title('simulation')
    
    xlabel('spike length (sec)')
    ylabel('count')
end

set(gcf,'Position',[100 100 1000 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrfxhrf_offset_baseline_sims.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% what if one stimulus onset is sustained

TL = tiledlayout(3,3,'TileIndexing','columnmajor') ; 

rng(42)

nn = 25 ; 

%randoffsets = (rand(nn,2)-0.5) .* 5 ;
randoffsets = normrnd(0,1,[nn 2]) ; 

hrf1 = getcanonicalhrf(1,0.72) ; 

for shiftn = 1:3

    nexttile([1 2])
    
    for idx = 1:nn
    
        hrf1 = getcanonicalhrf(rand(1).*7,0.72) ; 

        t1 = [zeros(shiftn,1) ; hrf1(1:60)' ] +randoffsets(idx,1) ; 

        hrf2 = getcanonicalhrf(rand(1).*7,0.72) ; 

        t2 = [ hrf2(1:60)' ; zeros(shiftn,1) ] +randoffsets(idx,2) ; 

        ee = prod([t1 t2],2) ;

        plot((0:length(ee)-1).*TR,ee,'Color',[ CM(20,:) 0.5],'LineWidth',2)
        hold on
        
    end
    hold off
    xlim([0 20])
    yline(2.25,'Color','b','LineWidth',1)
    
    title(['hrf x hrf with random start mag., offset ' num2str(shiftn.*TR) ' sec.'])

    if shiftn == 2
        xlabel('time (sec)')
    end

end

% lets do a simulation for each
nsim = 5000 ; 
spklen = zeros(nsim,1) ; 
%randoffsets = (rand(nsim,2)-0.5) .* 4 ; 
randoffsets = normrnd(0,1,[nsim 2]) ; 

histedges = [0:1:50].*TR ; 

for shiftn = 1:3

    nexttile
    
    for idx = 1:nsim
    
        hrf1 = getcanonicalhrf(rand(1).*7,0.72) ; 

        t1 = [zeros(shiftn,1) ; hrf1(1:60)' ] +randoffsets(idx,1) ; 

        hrf2 = getcanonicalhrf(rand(1).*7,0.72) ; 

        t2 = [ hrf2(1:60)' ; zeros(shiftn,1) ] +randoffsets(idx,2) ; 

        ee = prod([t1 t2],2) ;
        %ee = ee(10:30) ; % only focus on hrf rise
    
        [~,ss ] = spk_lenmat(ee>2.25) ; 
        ss=cell2mat(ss) ; 

        % if length(ss) > 1
        %     disp('feff')
        %     ss
        % end

        if ~isnan(ss)
            spklen(idx) =  max(ss) ; % max, in case theres more than 1 for some reason 
        end
    end
    
    % there are bound to be some weird trials where above 0 the whole
    % time... just keep them in, they don't affect much
    % spklen(snpkle>20) = nan ;

    spklen = spklen .* TR ;
    nt = sum(spklen>0) ;

    histogram(nonzeros(spklen(spklen<(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]) ; 
    hold on
    histogram(nonzeros(spklen(spklen>=(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]) ; 
    hold off

    xticks([5 10 15])
    xticklabels({ '5' '10' '15' })

    xlim([1 20])

    % comput hitrate >= 2.88 sec
    text(0.8,0.8,{ 'events >= 2.88' [ num2str(round(sum(spklen>=(4.*TR)) ./ nt,2)) '%']} , ...
        'units','normalized','HorizontalAlignment','center') ; 

    title('simulation')
    
    xlabel('spike length (sec)')
    ylabel('count')
end

set(gcf,'Position',[100 100 1000 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrfxhrf_offset_baseline_sustained_sims.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%

% ts = zscore(datStr(42).ts(:,1:200)) ; 
% ets = get_ets(ts) ; 
% 
% [~,ss1] = spk_lenmat(ts>prctile(ts,80)) ; 
% [~,ss2] = spk_lenmat(ets>prctile(ets,80)) ; 
% 
% histogram(cell2mat(ss1'),'Normalization','probability')
% hold on
% histogram(cell2mat(ss2'),'Normalization','probability')
% ylabel('probability')
% xlabel('spk length')
% 
% legend({'times series spk len' 'edge ts spk len'})
% 
% %%
% 
% rr = zeros(100,1) ;  
% for idx = 1:100
%     disp(idx)
% 
%     ts = zscore(datStr(idx).ts(:,1:200)) ; 
%     ets = get_ets(ts) ; 
% 
%     [~,ss1] = spk_lenmat(ts>prctile(ts,80)) ; 
%     [~,ss2] = spk_lenmat(ets>prctile(ets,80)) ; 
% 
%     rr(idx) = mean(cell2mat(ss1'))-mean(cell2mat(ss2')) ; 
% 
% end
% 
% histogram(rr)
% xlabel('node - edge hist. dist.')



