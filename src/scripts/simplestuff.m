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
spklen1 = zeros(nsim,1) ; 
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
            spklen1(idx) =  max(ss) ; % max, in case theres more than 1 for some reason 
        end
    end
    
    % there are bound to be some weird trials where above 0 the whole
    % time... just keep them in, they don't affect much
    % spklen(snpkle>20) = nan ;

    spklen1 = spklen1 .* TR ;
    nt = sum(spklen1>0) ;

    histogram(nonzeros(spklen1(spklen1<(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]) ; 
    hold on
    histogram(nonzeros(spklen1(spklen1>=(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]) ; 
    hold off

    xticks([5 10 15])
    xticklabels({ '5' '10' '15' })

    xlim([1 20])

    % comput hitrate >= 2.88 sec
    text(0.8,0.8,{ 'events >= 2.88' [ num2str(round(sum(spklen1>=(4.*TR)) ./ nt,2)) '%']} , ...
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
spklen1 = zeros(nsim,1) ; 
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
            spklen1(idx) =  max(ss) ; % max, in case theres more than 1 for some reason 
        end
    end
    
    % there are bound to be some weird trials where above 0 the whole
    % time... just keep them in, they don't affect much
    % spklen(snpkle>20) = nan ;

    spklen1 = spklen1 .* TR ;
    nt = sum(spklen1>0) ;

    histogram(nonzeros(spklen1(spklen1<(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.8 0.8 0.8]) ; 
    hold on
    histogram(nonzeros(spklen1(spklen1>=(4*TR))),histedges,...
        'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]) ; 
    hold off

    xticks([5 10 15])
    xticklabels({ '5' '10' '15' })

    xlim([1 20])

    % comput hitrate >= 2.88 sec
    text(0.8,0.8,{ 'events >= 2.88' [ num2str(round(sum(spklen1>=(4.*TR)) ./ nt,2)) '%']} , ...
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


%% do a grid test 

% hrflens = 1:10 ; 
% nlens = length(hrflens) ; 
% shiftlens = 0:5 ; 
% nshifts = length(shiftlens) ; 
% noiselevel = 0.25 ; 
% baselineshift = 0.25 ; 
% 
% nsim = 1000 ; 
% hrfamp = 2 ; 
% %wanttime = 1:30 ; 
% wanttime = find((getcanonicalhrf(1,0.72)*hrfamp).^2  > 2.25) ; 
% % pad the ends
% wanttime = [ wanttime(1)-1 wanttime wanttime(end)+1 ] ;
% 
% eventthr = 3 ; 
% 
% res1 = zeros(nlens,nshifts) ; 
% res2 = zeros(nlens,nshifts) ; 
% 
% for idx = 1:nlens
%     for jdx = 1:nshifts
% 
%         disp([num2str(idx) '-' num2str(jdx)])
%         shiftn = shiftlens(jdx) ; 
% 
%         spklen1 = zeros(nsim,1) ; 
%         spklen2 = zeros(nsim,1) ; 
% 
%         for sdx = 1:nsim
% 
%             hrf1 = getcanonicalhrf(1,TR).*hrfamp  ; 
%             hrf2 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 
% 
%             hrf1 = hrf1 + normrnd(0,noiselevel,size(hrf1)) + ((rand(1)-0.5)*baselineshift) ; 
%             hrf2 = hrf2 + normrnd(0,noiselevel,size(hrf2)) + ((rand(1)-0.5)*baselineshift) ; 
% 
%             % hrf1 = [zeros(1,shiftn) hrf1 ] ; 
%             % hrf2 = [ hrf2 zeros(1,shiftn) ] ; 
%             hrf1 = [ hrf1 zeros(1,shiftn) ] ; 
%             hrf2 = [ zeros(1,shiftn) hrf2  ] ; 
% 
%             nullrf  = generate_phase_surrogates(hrf2',0)' ; 
%             nullrf = (nullrf./(max(nullrf))).*hrfamp ; 
%             nullrf = nullrf + normrnd(0,noiselevel,size(nullrf)) + ((rand(1)-0.5)*baselineshift) ; 
% 
%             wt = wanttime ; 
% 
%             ee = prod([hrf1(wt)' hrf2(wt)'],2) ;
%             [~,ss ] = spk_lenmat(ee>eventthr) ; 
%             ss=cell2mat(ss) ; 
% 
%             if ~isnan(ss)
%                 spklen1(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
%             end
% 
%             ee = prod([hrf1(wt)' nullrf(wt)'],2) ;
%             [~,ss ] = spk_lenmat(ee>2.25) ; 
%             ss=cell2mat(ss) ; 
% 
%             if ~isnan(ss)
%                 spklen2(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
%             end
% 
%         end
% 
%         spklen1 = spklen1 .* TR ;
%         nt = sum(spklen1>0) ;
%         res1(idx,jdx) = sum(spklen1>=(4.*TR)) ./ nt ; 
% 
% 
%         spklen2 = spklen2 .* TR ;
%         nt = sum(spklen2>0) ;
%         res2(idx,jdx) = sum(spklen2>=(4.*TR)) ./ nt ; 
% 
%     end
% end
% 
% tiledlayout(1,2)
% 
% nexttile
% imagesc(res1)
% ylabel('hrf duration')
% xlabel('temporal shift')
% xticks(1:nshifts)
% xticklabels(cellstr(num2str(shiftlens')))
% cb = colorbar ; 
% cb.Label.String = 'hit rate' ; 
% clim([0 1])
% title('simulated event detection')
% 
% nexttile
% imagesc(res2)
% ylabel('hrf duration')
% xlabel('temporal shift')
% xticklabels(cellstr(num2str(shiftlens')))
% cb = colorbar ; 
% cb.Label.String = 'hit rate' ; 
% clim([0 1])
% title('simulated null event detection')

%% Both HRFs change their length

hrflens = 1:10 ; 
nlens = length(hrflens) ; 
shiftlens = 0.2:0.2:2 ; 
nshifts = length(shiftlens) ; 
noiselevels = [ 0.25 0.5 1] ;
nnoise = length(noiselevels) ; 
%baselineshift = 0.25 ; 

nsim = 500 ; 
hrfamp = 2 ; 
wanttime = 1:30 ; 
% wanttime = find((getcanonicalhrf(1,0.72)*hrfamp).^2  > 2.25) ; 
% % pad the ends
% wanttime = [ wanttime(1)-1 wanttime wanttime(end)+1 ] ;

eventthr = 2.25 ; 

res1 = zeros(nnoise,nlens,nshifts) ; 
res2 = zeros(nnoise,nlens,nshifts) ; 

rng(42)
for ndx = 1:nnoise

noiselevel = noiselevels(ndx) ; 

for idx = 1:nlens
    for jdx = 1:nshifts

        disp([num2str(idx) '-' num2str(jdx)])
        % shiftn = shiftlens(jdx) ; 

        spklen1 = zeros(nsim,1) ; 
        spklen2 = zeros(nsim,1) ; 

        for sdx = 1:nsim

            hrf1 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 
            hrf2 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 

            hrf2preshift = hrf2 + normrnd(0,noiselevel,size(hrf2)) ;
            hrf1 = hrf1 + normrnd(0,noiselevel,size(hrf1)) + ((rand(1)-0.5)*shiftlens(jdx)) ; 
            hrf2 = hrf2preshift + ((rand(1)-0.5)*shiftlens(jdx)) ; 

            nullrf  = generate_phase_surrogates(hrf2preshift',0)' ; 
            nullrf = (nullrf./(max(nullrf))).*hrfamp ; 
            nullrf = nullrf + normrnd(0,noiselevel,size(nullrf)) + ((rand(1)-0.5)*shiftlens(jdx)) ; 

            wt = wanttime ; 

            ee = prod([hrf1(wt)' hrf2(wt)'],2) ;
            [~,ss ] = spk_lenmat(ee>eventthr) ; 
            ss=cell2mat(ss) ; 

            if ~isnan(ss)
                spklen1(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
            end

            ee = prod([hrf1(wt)' nullrf(wt)'],2) ;
            [~,ss ] = spk_lenmat(ee>2.25) ; 
            ss=cell2mat(ss) ; 

            if ~isnan(ss)
                spklen2(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
            end

        end

        spklen1 = spklen1 .* TR ;
        nt = sum(spklen1>0) ;
        res1(ndx,idx,jdx) = sum(spklen1>=(4.*TR)) ./ nt ; 


        spklen2 = spklen2 .* TR ;
        nt = sum(spklen2>0) ;
        res2(ndx,idx,jdx) = sum(spklen2>=(4.*TR)) ./ nt ; 

    end
end

end

%%

tiledlayout(2,3,'TileIndexing','columnmajor')

for ndx = 1:3

nt=nexttile
h = imagesc(squeeze(res1(ndx,:,:))) ; 
ylabel('hrf duration')
xlabel('baseline var. mag')
xticks(1:nshifts)
xticklabels(cellstr(num2str(shiftlens')))
yticks(hrflens)
cb = colorbar ; 
cb.Label.String = 'hit rate' ; 
clim([0 1])
title(['simulated event detection, noise: ' num2str(noiselevels(ndx))])
colormap(nt,parula(100))

nt = nexttile

imagesc(squeeze(res2(ndx,:,:)))
ylabel('hrf duration')
xlabel('baseline var. mag')
xticks(1:nshifts)
xticklabels(cellstr(num2str(shiftlens')))
yticks(hrflens)
cb = colorbar ; 
cb.Label.String = 'hit rate' ; 
clim([0 1])
title(['simulated null event detection, noise: ' num2str(noiselevels(ndx))])

% nt = nexttile
% %calc F1 score 
% r1 = squeeze(res1(ndx,:,:)) ;
% r2 = squeeze(res2(ndx,:,:)) ; 
% 
% f1mat = (2.* (r1.*nsim)) ./  ...
%     ( (2.* (r1.*nsim)) + (r2.*nsim) + ((1-r1).*nsim) ) ; 
% 
% h = imagesc(f1mat) ; 
% ylabel('hrf duration')
% xlabel('baseline var. mag')
% xticks(1:nshifts)
% xticklabels(cellstr(num2str(shiftlens')))
% yticks(hrflens)
% cb = colorbar ; 
% cb.Label.String = 'f1' ; 
% clim([0 1])
% title('f1 score')
% colormap(nt,fblue(100))

end

%%

set(gcf,'Position',[100 100 1200 600])
set(gcf,'Color','w')
orient(gcf,'landscape')

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrfxhrf_grid_noise.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%% %% same as above but keep 1 hrf of length 1

hrflens = 1:10 ; 
nlens = length(hrflens) ; 
shiftlens = 0.2:0.2:2 ; 
nshifts = length(shiftlens) ; 
noiselevels = [ 0.25 0.5 1] ;
nnoise = length(noiselevels) ; 
%baselineshift = 0.25 ; 

nsim = 500 ; 
hrfamp = 2 ; 
wanttime = 1:30 ; 
% wanttime = find((getcanonicalhrf(1,0.72)*hrfamp).^2  > 2.25) ; 
% % pad the ends
% wanttime = [ wanttime(1)-1 wanttime wanttime(end)+1 ] ;

eventthr = 2.25 ; 

res1 = zeros(nnoise,nlens,nshifts) ; 
res2 = zeros(nnoise,nlens,nshifts) ; 

rng(42)
for ndx = 1:nnoise

noiselevel = noiselevels(ndx) ; 

for idx = 1:nlens
    for jdx = 1:nshifts

        disp([num2str(idx) '-' num2str(jdx)])
        % shiftn = shiftlens(jdx) ; 

        spklen1 = zeros(nsim,1) ; 
        spklen2 = zeros(nsim,1) ; 

        for sdx = 1:nsim

            hrf1 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 
            hrf2 = [ zeros(1,hrflens(idx)) getcanonicalhrf(1,TR).*hrfamp ]   ; 

            hrf2preshift = hrf2 + normrnd(0,noiselevel,size(hrf2)) ;
            hrf1 = hrf1 + normrnd(0,noiselevel,size(hrf1)) + ((rand(1)-0.5)*shiftlens(jdx)) ; 
            hrf2 = hrf2preshift + ((rand(1)-0.5)*shiftlens(jdx)) ; 

            nullrf  = generate_phase_surrogates(hrf2preshift',0)' ; 
            nullrf = (nullrf./(max(nullrf))).*hrfamp ; 
            nullrf = nullrf + normrnd(0,noiselevel,size(nullrf)) + ((rand(1)-0.5)*shiftlens(jdx)) ; 

            wt = wanttime ; 

            ee = prod([hrf1(wt)' hrf2(wt)'],2) ;
            [~,ss ] = spk_lenmat(ee>eventthr) ; 
            ss=cell2mat(ss) ; 

            if ~isnan(ss)
                spklen1(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
            end

            ee = prod([hrf1(wt)' nullrf(wt)'],2) ;
            [~,ss ] = spk_lenmat(ee>2.25) ; 
            ss=cell2mat(ss) ; 

            if ~isnan(ss)
                spklen2(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
            end

        end

        spklen1 = spklen1 .* TR ;
        nt = sum(spklen1>0) ;
        res1(ndx,idx,jdx) = sum(spklen1>=(4.*TR)) ./ nt ; 


        spklen2 = spklen2 .* TR ;
        nt = sum(spklen2>0) ;
        res2(ndx,idx,jdx) = sum(spklen2>=(4.*TR)) ./ nt ; 

    end
end

end

%%

tiledlayout(2,3,'TileIndexing','columnmajor')

for ndx = 1:3

nt=nexttile
h = imagesc(squeeze(res1(ndx,:,:))) ; 
ylabel('hrf duration')
xlabel('baseline var. mag')
xticks(1:nshifts)
xticklabels(cellstr(num2str(shiftlens')))
yticks(hrflens)
cb = colorbar ; 
cb.Label.String = 'hit rate' ; 
clim([0 1])
title(['simulated event detection, noise: ' num2str(noiselevels(ndx))])
colormap(nt,parula(100))

nt = nexttile

% % calc F1 score 
% r1 = squeeze(res1(ndx,:,:)) ;
% r2 = squeeze(res2(ndx,:,:)) ; 
% 
% f1mat = (2.* (r1.*nsim)) ./  ...
%     ( (2.* (r1.*nsim)) + (r2.*nsim) + ((1-r1).*nsim) ) ; 
% 
% h = imagesc(f1mat) ; 
% ylabel('hrf duration')
% xlabel('baseline var. mag')
% xticks(1:nshifts)
% xticklabels(cellstr(num2str(shiftlens')))
% yticks(hrflens)
% cb = colorbar ; 
% cb.Label.String = 'f1' ; 
% clim([0 1])
% title('f1 score')
% colormap(nt,fblue(100))

imagesc(squeeze(res2(ndx,:,:)))
ylabel('hrf duration')
xlabel('baseline var. mag')
xticks(1:nshifts)
xticklabels(cellstr(num2str(shiftlens')))
yticks(hrflens)
cb = colorbar ; 
cb.Label.String = 'hit rate' ; 
clim([0 1])
title(['simulated null event detection, noise: ' num2str(noiselevels(ndx))])

end

set(gcf,'Position',[100 100 1200 600])
set(gcf,'Color','w')

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/hrfxhrf_grid_noise.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)


%% DOESNT WORK ADDING THE NOISE THIS WAY
% 
% tiledlayout(2,3,'TileIndexing','columnmajor')
% 
% hrflens = 1:10 ; 
% nlens = length(hrflens) ; 
% shiftlens = 0.2:0.2:2 ; 
% nshifts = length(shiftlens) ; 
% %noiselevel = 0.25  ;
% %baselineshift = 0.25 ; 
% 
% nsim = 500 ; 
% hrfamp = 2 ; 
% wanttime = 1:30 ; 
% % wanttime = find((getcanonicalhrf(1,0.72)*hrfamp).^2  > 2.25) ; 
% % % pad the ends
% % wanttime = [ wanttime(1)-1 wanttime wanttime(end)+1 ] ;
% 
% noiselevels = [ 0.5 1 1.5 ] ; 
% 
% % eventthrs = [1.75 2.25 2.75] ; 
% % nthr = length(eventthrs) ; 
% eventthr = 2.25 ;
% 
% res1 = zeros(nnoise,nlens,nshifts) ; 
% res2 = zeros(nnoise,nlens,nshifts) ; 
% 
% for ndx = 1:nthr
% 
% noiselevel = noiselevels(ndx) ; 
% 
% for idx = 1:nlens
%     for jdx = 1:nshifts
% 
%         disp([num2str(idx) '-' num2str(jdx)])
%         % shiftn = shiftlens(jdx) ; 
% 
%         spklen1 = zeros(nsim,1) ; 
%         spklen2 = zeros(nsim,1) ; 
% 
%         for sdx = 1:nsim
% 
%             hrf1 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 
%             hrf2 = getcanonicalhrf(hrflens(idx),TR).*hrfamp  ; 
% 
%             n1 = generate_phase_surrogates(hrf1',0)' ;
%             n2 = generate_phase_surrogates(hrf2',0)' ;
% 
%             n1 = (n1./max(abs(n1))).*noiselevel; 
%             n2 = (n2./max(abs(n2))).*noiselevel; 
% 
%             hrf1 = hrf1 + n1 + ((rand(1)-0.5)*shiftlens(jdx)) ; 
%             hrf2 = hrf2 + n2 + ((rand(1)-0.5)*shiftlens(jdx)) ; 
% 
%             % hrf2preshift = hrf2 + normrnd(0,noiselevel,size(hrf2)) ;
%             % hrf1 = hrf1 + normrnd(0,noiselevel,size(hrf1)) + ((rand(1)-0.5)*shiftlens(jdx)) ; 
%             % hrf2 = hrf2preshift + ((rand(1)-0.5)*shiftlens(jdx)) ; 
% 
%             nullrf  = generate_phase_surrogates(getcanonicalhrf(hrflens(idx),TR)')' ; 
%             nullrf = (nullrf./(max(nullrf))).*(rand(1).*hrfamp) ; 
%             nullrf = nullrf + ((rand(1)-0.5)*shiftlens(jdx)) ; 
%             % nullrf = -hrf2 ; 
% 
%             wt = wanttime ; 
% 
%             ee = prod([hrf1(wt)' hrf2(wt)'],2) ;
%             [~,ss ] = spk_lenmat(ee>eventthr) ; 
%             ss=cell2mat(ss) ; 
% 
%             if ~isnan(ss)
%                 spklen1(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
%             end
% 
%             ee = prod([hrf1(wt)' nullrf(wt)'],2) ;
%             [~,ss ] = spk_lenmat(ee>2.25) ; 
%             ss=cell2mat(ss) ; 
% 
%             if ~isnan(ss)
%                 spklen2(sdx) =  max(ss) ; % max, in case theres more than 1 for some reason 
%             end
% 
%         end
% 
%         spklen1 = spklen1 .* TR ;
%         nt = sum(spklen1>0) ;
%         res1(ndx,idx,jdx) = sum(spklen1>=(4.*TR)) ./ nt ; 
% 
% 
%         spklen2 = spklen2 .* TR ;
%         nt = sum(spklen2>0) ;
%         res2(ndx,idx,jdx) = sum(spklen2>=(4.*TR)) ./ nt ; 
% 
%     end
% end
% 
% 
% nexttile
% imagesc(squeeze(res1(ndx,:,:)))
% ylabel('hrf duration')
% xlabel('baseline var. mag')
% xticks(1:nshifts)
% xticklabels(cellstr(num2str(shiftlens')))
% yticks(hrflens)
% cb = colorbar ; 
% cb.Label.String = 'hit rate' ; 
% clim([0 1])
% title('simulated event detection')
% 
% nexttile
% 
% % calc F1 score 
% r1 = squeeze(res1(ndx,:,:)) ;
% r2 = squeeze(res2(ndx,:,:)) ; 
% 
% f1mat = (2.* (r1.*nsim)) ./  ...
%     ( (2.* (r1.*nsim)) + (r2.*nsim) + ((1-r1).*nsim) ) ; 
% 
% 
% imagesc(f1mat)
% ylabel('hrf duration')
% xlabel('baseline var. mag')
% xticks(1:nshifts)
% xticklabels(cellstr(num2str(shiftlens')))
% yticks(hrflens)
% cb = colorbar ; 
% cb.Label.String = 'f1 score' ; 
% clim([0 1])
% title('simulated null event detection')
% 
% end

%%

hrf = getcanonicalhrf(TR,1) ; 

bc = zeros(100,1) ; 
bc(2:10) = 1 ; 

[~,cc] = pad_conv(bc,hrf,100) ; 

%%

writematrix([0.1:0.2:20]','/Users/faskowitzji/joshstuff/sandbox/hrfmag/evalhrftimes.txt')

mat = zeros(100,100) ; 
for idx = 1:100
    mat(:,idx) = load([ '/Users/faskowitzji/joshstuff/sandbox/hrfmag/len' num2str(idx) '.txt' ] ) ; 
end
