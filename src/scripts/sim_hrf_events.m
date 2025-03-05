%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_read10.m') 

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

%%

tiledlayout(2,1)

spklens = 1:2:29 ; 
CM = internet(max(spklens)) ; 

nexttile

for ii = spklens

    hh = getcanonicalhrf_wmag(ii,TR) ; 
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

ylabel('amplitude')

%cb = colorbar() ; 

nexttile
for ii = spklens

    % hh = getcanonicalhrf(ii,TR) ; 
    % ll = length(hh) ; 
    % % do a fair zscore, relfecting baseline at '0'
    % tmp = zscore([ hh -hh]) ; 
    % hh = tmp(1:ll) ; 
    % 
    % hhh = prod([hh' hh'],2) ;
    % 
    % plot((0:length(hhh)-1).*0.72,hhh,'Color',CM(ii,:),'LineWidth',2)
    % xlim([0 45])
    % hold on

    hh = getcanonicalhrf_wmag(ii,TR) ; 

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


hrf1 = getcanonicalhrf_wmag(1,TR) ; 
ll = length(hrf1) ; 
% do a fair zscore, relfecting baseline at '0'
% tmp = zscore([ hrf1 -hrf1]) ; 
% hrf1 = tmp(1:ll) ; 

nexttile()

plot((0:ll-1).*TR,getcanonicalhrf_wmag(1,TR),'LineWidth',2)
xlim([0 30])

title('hrf 1 sec. stim.')

nexttile()

% plot((0:ll-1).*TR,prod([ hrf1' hrf1' ],2),'LineWidth',2)

for idx = 0:10

    %hh1 = [  hrf1'  hrf1' ] ; 
    hh2 = [ [ hrf1' ; zeros(idx,1)] [ zeros(idx,1) ; hrf1'] ] ; 
    % tmp = zscore( [hh2 ; -hh2 ]) ; 
    % hh2 = tmp(1:ll,:) ; 
    

    hold on
    plot((0:ll-1).*TR,prod(hh2(1:ll,:),2),'LineWidth',2,'Color',CM(idx+5,:))
    
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

%%

TR = finfo.TR ; 
tiledlayout(1,3,"TileSpacing",'none') ;

hrfthrs = [2 2.25 2.5 ] ; 

for tdx = 1:3

    hrfthr = hrfthrs(tdx) ; 

    hrf = getcanonicalhrf_wmag(1,1) ; 
    
    hrf2_onsets = 0:0.1:3 ; 
    hrf2_lengths = 1:0.5:15 ; 
    % magDelta = 1:0.1:4 ; 
    
    res1 = zeros(length(hrf2_onsets),length(hrf2_lengths)) ; 
    
    for idx = 1:length(hrf2_onsets)
        for jdx = 1:length(hrf2_lengths)
    
            hrf1 = hrf.*1 ; 
            hrf2 = getcanonicalhrf_wmag(hrf2_lengths(jdx),1) ; 
    
            h1 = [ hrf1 ] ;
            
            h2 = interp1((0:length(hrf2)-1)+hrf2_onsets(idx),hrf2,(0:length(hrf2)-1),"linear") ; 
            h2 = h2(1:length(h1)) ; 
            h2(isnan(h2)) = 0 ; 
    
            ee = prod([h1(:) h2(:)],2) ; 
    
            res1(idx,jdx) = meas_abv_thr_interp(ee,hrfthr)*TR ; 
    
        end
    end

    nexttile()

    h = imagesc(res1)
    h.AlphaData = ~isnan(res1)
    title(['event thr: ' num2str(hrfthrs(tdx)) ' co-fluc.' ])
    clim([0 4])
    if tdx == 3
    cb = colorbar()
    cb.Label.String = 'event length (sec)'
    end

    if tdx == 1
    yticks(1:length(hrf2_onsets))
    yticklabels(num2str(hrf2_onsets'))
    ylabel('delta onset delta (sec)')
    else
        yticks([])
    end

    xticks(1:length(hrf2_lengths))
    xticklabels(num2str(hrf2_lengths'))
    xlabel('stim duration (sec)')

end

colormap(parula(100))

%% print it

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[0 100 1200 400])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/grid_deltaonset_n_len.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%

neuralHz = 10 ; 

poissonSpkGen = @(firingrateHz,lenSec) rand(1,floor(lenSec*neuralHz)) < (firingrateHz*1/neuralHz) ; 
regularSpkGen = @(firingrateHz,lenSec) ismember(1:lenSec*neuralHz,floor(linspace(1,lenSec*neuralHz,floor(firingrateHz*lenSec)))) ; 
irregularSpkGen = @(firingrateHz,lenSec,poissShiftVar) ...
    ismember(1:lenSec*neuralHz, arrayfun(@(x_) x_ + (poissrnd(poissShiftVar)*randsample([-1,1],1)) , find(regularSpkGen(firingrateHz,lenSec)) )) ; 
boxCarSpkGen = @(firingrateHz,lenSec) ones(floor(lenSec*neuralHz),1) ; 

hrf = getcanonicalhrf_wmag(1/neuralHz,1/neuralHz) ; 

useGenFunc = boxCarSpkGen ; 

%%

% addnoisefrac = 0.5 ; 
% fmriHz = 1/TR ; 
% 
% spkThrSweep = ( 2:0.05:3 ) ; 
% stimLenSweep = ( 1:0.5:10 ) ; 
% 
% nrep=500 ; 
% res1 = zeros(length(spkThrSweep),length(stimLenSweep),nrep) ; 
% res2 = zeros(length(spkThrSweep),length(stimLenSweep),nrep) ; 
% 
% rng(42)
% for ndx = 1:nrep
% 
% disp(ndx)
% for idx = 1:length(spkThrSweep)
%     for jdx = 1:length(stimLenSweep)
% 
%         %% neural ts get noise before convolution
% 
%         % init ts of 60 sec
%         nts1 = zeros(neuralHz*60,1) ; 
%         nts2 = zeros(neuralHz*60,1) ; 
% 
%         sig1Len = 1 ; 
% 
%         % insert the 'neural activity' at 1 sec 
%         insertInds1 = (2*neuralHz):((2*neuralHz)+(sig1Len*neuralHz)-1) ;
%         insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%         nts1(insertInds1) = useGenFunc(NaN,sig1Len) ; 
%         nts2(insertInds2) = useGenFunc(NaN,stimLenSweep(jdx)) ; 
% 
%         % convolve with the hrf at neuralHz
%         [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%         [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%         % resample to fmriHZ
%         fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%         fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%         % adding noise
%         % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%         % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%         % add 2% noise ~ assumption 50% of max signal amp
%         fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*addnoisefrac)) ; 
%         fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*addnoisefrac)) ; 
% 
%         ee = prod([fts1(:) fts2(:)],2) ; 
% 
%         mm = meas_abv_thr_interp(ee,spkThrSweep(idx)).*TR ; 
%         if mm<0
%             error('neg length problem')
%         end
% 
%         res2(idx,jdx,ndx) =  mm(1) ;
% 
%     end
% end
% 
% end
% 
% %
% 
% mat = mean(res2,3,'omitnan') ; 
% 
% % set noise to %2 of signal
% % do boxcar
% % take a few thresholds 
% 
% exemplarindz = [6 6; 6 10; 11 7 ] ; 
% exemplarindxcol = [0 0 1 ; 1 0 0 ; 0 1 0 ] ; 
% 
% TL = tiledlayout(1,1+size(exemplarindxcol,1)) ; 
% 
% nexttile(TL)
% 
% 
% %res2(isnan(res2)) = 0 ; 
% h = imagesc(mat) ; 
% h.AlphaData = ~isnan(mat) ; 
% cb = colorbar() ; 
% cb.Label.String = 'edge above thr. duration' ;
% 
% yticks(1:length(spkThrSweep))
% yticklabels(num2str(spkThrSweep'))
% ylabel('co-fluct. threshold')
% 
% xticks(1:length(stimLenSweep))
% xticklabels(num2str(stimLenSweep'))
% xlabel('neural activity duration')
% 
% axis square
% colormap(acton(100))
% 
% clim([0 max(mat,[],'all')])
% 
% for edx = 1:size(exemplarindxcol,1)
%     imsc_addsquare(1:length(spkThrSweep)==exemplarindz(edx,1), ...
%         1:length(stimLenSweep)== exemplarindz(edx,2), ...
%         0,exemplarindxcol(edx,:))
% end
% 
% for edx = 1:size(exemplarindxcol,1)
% 
% % get an exemplar
% 
% nreps = 50 ; 
% exres = cell(nreps,1) ; 
% 
% idx = exemplarindz(edx,1) ; 
% jdx = exemplarindz(edx,2) ; 
% 
% for ndx = 1:nreps
%     % init ts of 60 sec
%     nts1 = zeros(neuralHz*60,1) ; 
%     nts2 = zeros(neuralHz*60,1) ; 
% 
%     % insert the 'neural activity' at 1 sec 
%     sig1Len = 1 ; 
% 
%     insertInds1 = (2*neuralHz):((2*neuralHz)+(sig1Len*neuralHz)-1) ;
%     insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%     nts1(insertInds1) = useGenFunc(NaN,sig1Len) ; 
%     nts2(insertInds2) = useGenFunc(NaN,stimLenSweep(jdx)) ; 
% 
%     % convolve with the hrf at neuralHz
%     [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%     [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%     % resample to fmriHZ
%     fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%     fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%     % adding noise
%     % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%     % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%     % add 2% signal noise 
%     % fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*addnoisefrac)) ; 
%     % fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*addnoisefrac)) ; 
%     fts1 = fts1 + generate_phase_surrogates(fts1).*(max(fts1).*addnoisefrac) ; 
%     fts2 = fts2 + generate_phase_surrogates(fts2).*(max(fts2).*addnoisefrac) ; 
% 
%     ee = prod([fts1(:) fts2(:)],2) ; 
% 
%     exres{ndx} = [fts1 fts2 ee] ; 
% 
% end
% 
% %
% 
% TL2 = tiledlayout(TL,3,1) ; 
% TL2.Layout.Tile = edx+1 ;
% 
% nexttile(TL2,[2 1])
% 
% timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
% wt = 1:find(timevec>45,1,'first') ; 
% 
% for ndx = 1:nreps
% 
%     plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
%     hold on
%     plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])
% 
% 
%     plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])
% 
% 
% end
% 
% text(0.5,0.4,{ ['threshold: ' num2str(spkThrSweep(idx))] ...
%     ['neural duration: ' num2str(stimLenSweep(jdx)) ' sec.'] ...
%     ['abv thr dur.: ' num2str(mat(idx,jdx)) ' sec.']},'units','normalized')
% 
% xlim([0 45])
% ylim([-2 10])
% 
% title('example time series')
% legend({'node time series' '' 'edge time series'})
% 
% %
% 
% rng(42)
% 
% nex = 1 ; 
% exnts = zeros(neuralHz*60,nex) ; 
% 
% for tdx = 1:nex
% 
% % insert the 'neural activity' at 1 sec 
% insertInds = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
% exnts(insertInds,tdx) = useGenFunc(NaN,stimLenSweep(jdx)) ; 
% 
% end
% 
% nexttile(TL2,[1 1])
% 
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts+0.5),'Color','b')
% % plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
% ylim([0 2])
% xlim([0 45])
% 
% yticks([])
% 
% title('example boxcar')
% xlabel('time (sec)')
% 
% end
% 
% %% print it
% 
% set(gcf,'Color','w')
% orient(gcf,'landscape')
% set(gcf,'Position',[0 100 1500 300])
% 
% out_figdir = [ './reports/figures/supp/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/grid_thr_by_len.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)


%% now vary the length of both durations 

addvarfrac=0.5 ; 
%addnoisefrac = 0.4 ; 
snrtarg = 2 ; 

hrf = getcanonicalhrf_wmag(1/neuralHz,1/neuralHz) ; 

fmriHz = 1/finfo.TR ; 

%spkThrSweep = ( 2:0.05:3 ) ; 
stimLenSweep = ( 0.1:0.2:5 ) ; 

nrep=500 ; 
res2 = zeros(length(stimLenSweep),length(stimLenSweep),nrep) ; 
res2null = zeros(length(stimLenSweep),length(stimLenSweep),nrep) ; 
cnr2 = zeros(length(stimLenSweep),length(stimLenSweep),nrep) ;

%%

rng(42)
for ndx = 1:nrep

disp(ndx)
for idx = 1:length(stimLenSweep)
    for jdx = 1:length(stimLenSweep)

        %% neural ts get noise before convolution

        % init ts of 60 sec
        nts1 = zeros(neuralHz*60,1) ; 
        nts2 = zeros(neuralHz*60,1) ; 

        %sig1Len = 1 ; 

        % insert the 'neural activity' at 1 sec 
        insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
        insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;

        nts1(insertInds1) = 1 ; 
        nts2(insertInds2) = 1 ; 

        % convolve with the hrf at neuralHz
        [~,cnts1] = pad_conv(nts1,hrf,100) ; 
        [~,cnts2] = pad_conv(nts2,hrf,100) ; 

        % resample to fmriHZ
        
        fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
        fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
        
        % % adding noise
        % % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
        % % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
        % 
        % % add variability 
        % % fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*addnoisefrac)) ; 
        % % fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*addnoisefrac)) ; 
        % nn1 = generate_phase_surrogates(fts1); 
        % nn2 = generate_phase_surrogates(fts2) ; 
        % nn1 = (nn1./max(nn1)).*(max(fts1).*addvarfrac) ; 
        % nn2 = (nn2./max(nn2)).*(max(fts2).*addvarfrac) ; 
        % 
        % cnr2(idx,jdx,ndx) = mean([ max(fts1)/std(nn1) max(fts2)/std(nn2)]) ; 
        % 
        % fts1 = fts1 + nn1 ; 
        % fts2 = fts2 + nn2 ; 

        fts1 = add_var_n_noiseSNR(fts1,addvarfrac,snrtarg) ; 
        fts2 = add_var_n_noiseSNR(fts2,addvarfrac,snrtarg) ; 

        ee = prod([fts1(:) fts2(:)],2) ; 

        mm = meas_abv_thr_interp(ee,2.25).*TR ; 
        if mm<0
            error('neg length problem')
        end

        res2(idx,jdx,ndx) =  mm(1) ;

        %% now the null

        nfts1 = generate_phase_surrogates(fts1) ; 
        nfts2 = generate_phase_surrogates(fts2) ; 

        nee = prod([nfts1(:) nfts2(:)],2) ; 

        % 1-21 means an event during the first 15 sec
        mm = meas_abv_thr_interp(nee(1:21),2.25).*TR ; 
        mm = mm(1) ; 
        if mm<0
            error('neg length problem')
        end

        res2null(idx,jdx,ndx) = mm(1) ; 

    end
end

end

durrmat = mean(res2,3,'omitnan') ; 
nullmat = mean(res2null,3,'omitnan') ; 

%% save the info

filename = [ DD.PROC '/sim_hrf_grid_' OUTSTR '.mat' ] ; 
save(filename,'res2','res2null','durrmat','nullmat','-v7.3')

%% VIZ part

rng(42)

exemplarindz = [23 6; 15 16; 11 16 ] ; 
exemplarindxcol = [0 0 1 ; 1 0 0 ; 0 1 0 ] ; 

TL = tiledlayout(2,1+size(exemplarindxcol,1),'TileIndexing','rowmajor') ; 

nexttile(TL)

h = imagesc(durrmat) ; 
h.AlphaData = ~isnan(durrmat) ; 
cb = colorbar() ; 
cb.Label.String = 'edge above thr. duration' ;

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')

axis square
colormap(acton(100))

clim([0 max(durrmat,[],'all')])

% high entries
squrents = ~isnan(durrmat) & (durrmat >= 2.88 ) & (durrmat < 5.04) ; 
[u,v] = find(squrents) ; 
n = size(durrmat,1) ; 
for edx = 1:length(u)

    % imsc_addsquare(ca1,ca2,mirrorit,linecolor,linewidth)

    imsc_addsquare(1:n==u(edx), ...
        1:n==v(edx), ...
        0,[ 1 1 1 ],2)

end

% plot exemplars
for edx = 1:size(exemplarindxcol,1)
    imsc_addsquare(1:length(stimLenSweep)==exemplarindz(edx,1), ...
        1:length(stimLenSweep)== exemplarindz(edx,2), ...
        0,exemplarindxcol(edx,:))
end

for edx = 1:size(exemplarindxcol,1)

% get an exemplar

nreps = 20 ; 
exres = cell(nreps,1) ; 

idx = exemplarindz(edx,1) ; 
jdx = exemplarindz(edx,2) ; 

for ndx = 1:nreps
    % init ts of 60 sec
    nts1 = zeros(neuralHz*60,1) ; 
    nts2 = zeros(neuralHz*60,1) ; 

    % insert the 'neural activity' at 1 sec 
    sig1Len = 1 ; 

    insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
    insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;

    nts1(insertInds1) = 1 ; 
    nts2(insertInds2) = 1 ; 

    % convolve with the hrf at neuralHz
    [~,cnts1] = pad_conv(nts1,hrf,100) ; 
    [~,cnts2] = pad_conv(nts2,hrf,100) ; 

    % resample to fmriHZ
    fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
    fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 

    % adding noise
    fts1 = add_var_n_noiseSNR(fts1,addvarfrac,snrtarg) ; 
    fts2 = add_var_n_noiseSNR(fts2,addvarfrac,snrtarg) ; 


    ee = prod([fts1(:) fts2(:)],2) ; 

    exres{ndx} = [fts1 fts2 ee] ; 

end

TL2 = tiledlayout(TL,3,1) ; 
TL2.Layout.Tile = edx+1 ;

nexttile(TL2,[2 1])

timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
wt = 1:find(timevec>45,1,'first') ; 

for ndx = 1:nreps

    plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
    hold on
    plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])


    plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])


end

text(0.5,0.4,{ ...
    ['duration A: '  num2str(stimLenSweep(idx)) ' sec' ] ...
    ['duration B: ' num2str(stimLenSweep(jdx)) ' sec'] ...
    ['abv. thr. dur.: ' num2str(durrmat(idx,jdx)) ' sec'] }, ...
    'units','normalized')


xlim([0 45])
ylim([-2 10])

title('example time series')
legend({'node time series' '' 'edge time series'})

%nex = 1 ; 
exnts = zeros(neuralHz*60,2) ; 

% insert the 'neural activity' at 1 sec 
insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;

exnts(insertInds1,1) = 1 ; 
exnts(insertInds2,2) = 1 ; 


nexttile(TL2,[1 1])

plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+fliplr((0.5:2)),'Color','b')
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
ylim([0 2.5])
xlim([0 45])

yticks([.5 1.5])
yticklabels({'B' 'A'})

title('example boxcar')
xlabel('time (sec)')

end

%

nexttile(TL)

%res2(isnan(res2)) = 0 ; 

h = imagesc(nullmat) ; 
h.AlphaData = ~isnan(nullmat) ; 
cb = colorbar() ; 
cb.Label.String = 'edge above thr. duration' ;

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')

axis square

clim([0 max(durrmat,[],'all')])

for edx = 1:size(exemplarindxcol,1)
    imsc_addsquare(1:length(stimLenSweep)==exemplarindz(edx,1), ...
        1:length(stimLenSweep)== exemplarindz(edx,2), ...
        0,exemplarindxcol(edx,:))
end

for edx = 1:size(exemplarindxcol,1)

% get an exemplar

nreps = 20 ; 
exres = cell(nreps,1) ; 

idx = exemplarindz(edx,1) ; 
jdx = exemplarindz(edx,2) ; 

for ndx = 1:nreps
    % init ts of 60 sec
    nts1 = zeros(neuralHz*60,1) ; 
    nts2 = zeros(neuralHz*60,1) ; 

    % insert the 'neural activity' at 1 sec 
    sig1Len = 1 ; 

    insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
    insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;

    nts1(insertInds1) = 1 ; 
    nts2(insertInds2) = 1 ; 

    % convolve with the hrf at neuralHz
    [~,cnts1] = pad_conv(nts1,hrf,100) ; 
    [~,cnts2] = pad_conv(nts2,hrf,100) ; 

    % resample to fmriHZ
    fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
    fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 

    % adding noise
    % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
    % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;

    % add 2% signal noise 
    % nn1 = generate_phase_surrogates(fts1) ; 
    % nn2 = generate_phase_surrogates(fts2); 
    % nn1 = (nn1./max(nn1)).*(max(fts1).*addvarfrac) ; 
    % nn2 = (nn2./max(nn2)).*(max(fts2).*addvarfrac) ; 
    % 
    % fts1 = fts1 + nn1 ; 
    % fts2 = fts2 + nn2 ; 


    fts1 = add_var_n_noiseSNR(fts1,addvarfrac,snrtarg) ; 
    fts2 = add_var_n_noiseSNR(fts2,addvarfrac,snrtarg) ; 

    nfts1 = generate_phase_surrogates(fts1);
    nfts2 = generate_phase_surrogates(fts2);
 
    nee = prod([nfts1(:) nfts2(:)],2) ; 

    exres{ndx} = [nfts1 nfts2 nee] ; 

end

TL2 = tiledlayout(TL,3,1) ; 
TL2.Layout.Tile = edx+5 ;

nexttile(TL2,[2 1])

timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
wt = 1:find(timevec>45,1,'first') ; 

for ndx = 1:nreps

    plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
    hold on
    plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])


    plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])


end

hr = (nrep-sum(isnan(res2null),3))./nrep ; 

text(0.1,0.7,{ ...
    ['duration A: '  num2str(stimLenSweep(idx)) ' sec' ] ...
    ['duration B: ' num2str(stimLenSweep(jdx)) ' sec'] ...
    ['abv. thr. dur.: ' num2str(nullmat(idx,jdx)) ' sec'] ...
    ['hit rate: ' num2str(hr(idx,jdx)*100) '%' ] }, ...
    'units','normalized')


xlim([0 45])
ylim([-2 10])

title('example time series')
legend({'node time series' '' 'edge time series'})

%

rng(42)

%nex = 1 ; 
exnts = zeros(neuralHz*60,2) ; 

% insert the 'neural activity' at 1 sec 
insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;

exnts(insertInds1,1) = 1 ; 
exnts(insertInds2,2) = 1 ; 


nexttile(TL2,[1 1])

plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+fliplr((0.5:2)),'Color','b')
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
ylim([0 2.5])
xlim([0 45])

yticks([.5 1.5])
yticklabels({'B' 'A'})

title('example boxcar')
xlabel('time (sec)')

end

% mean(mean(cnr2,3),'all')
% 
% ans =
% 
%     4.5359

% std(mean(cnr2,3),[],'all')
% 
% ans =
% 
%     0.0312

%% print it

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[0 100 1500 600])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/grid_len_by_len.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% plot the HR matrix

hr = (nrep-sum(isnan(res2null),3))./nrep ; 

TL = tiledlayout(1,4)
TL.Title.String = 'Null data detection rates' ; 

nexttile
h = imagesc(hr.*100) ; 
h.AlphaData = hr~=0 ; 
cb = colorbar() ; 
cb.Label.String = 'above threshold event detection rate' ;
clim([0 100])
colormap(acton(100))
axis square

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')

% nexttile
% h = imagesc(hr.*nmat) ; 
% h.AlphaData = hr~=0 ; 
% cb = colorbar() ; 
% cb.Label.String = 'detection rate * mean abv. thr time' ;
% axis square
% 
% yticks(1:length(stimLenSweep))
% yticklabels(num2str(stimLenSweep'))
% ylabel('neural activity duration node A')
% 
% xticks(1:length(stimLenSweep))
% xticklabels(num2str(stimLenSweep'))
% xlabel('neural activity duration node B')

lll = stimLenSweep+stimLenSweep' ; 

cm = gray(100); 

nt = nexttile()
scatter(durrmat(:),nullmat(:),30,lll(:),'filled')
xlabel('above thr. duration, sim')
ylabel('above thr. duration, null')
xlim([0 9])
ylim([0 9])
refline(1,0)
colormap(nt,cm(10:90,:))
cb = colorbar() ; 
cb.Label.String = 'addivite activity dur.'
axis square
nt = nexttile



% scatter of time in simulated, vs hitrate for null
resmat = durrmat .* 1;
resmat(isnan(resmat)) = 0 ; 
h = scatter(resmat(:),hr(:).*100,30,lll(:),'filled') ; 
%h.MarkerFaceColor = [0.8 0.8 0.8 ] ; 
%h.MarkerFaceAlpha = 0.8 ;
colormap(nt,cm(10:90,:))
cb = colorbar() ; 
cb.Label.String = 'addivite activity dur.'
ylim([0 100])
xlim([0 10])

axis square
xlabel('simulated event length (non-null)')
ylabel('detection rate (null)')

nt = nexttile()

% scatter of time in simulated, vs hitrate for null
resmat = durrmat .* 1;
resmat(isnan(resmat)) = 0 ; 
h = scatter(resmat(resmat<5.04),hr(resmat<5.04).*100,30,lll(resmat<5.04),'filled') ; 
%h.MarkerFaceColor = [0.8 0.8 0.8 ] ; 
%h.MarkerFaceAlpha = 0.8 ;
colormap(nt,cm(10:90,:))
cb = colorbar() ; 
cb.Label.String = 'addivite activity dur.'
ylim([0 100])
axis square
xlabel('simulated event length (non-null)')
ylabel('detection rate (null)')
ylim([0 100])
xlim([0 10])
colormap(nt,cm(10:90,:))

yline(max(hr(resmat<5.04).*100))
text(0.6,max(hr(resmat<5.04).*100)*1.4/100,[ 'max: ' num2str(round(max(hr(resmat<5.04).*100),2)) ],'units','normalized')

%%

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[0 100 1000 400])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/grid_null_mats.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% make supp sim matrix with colored tiles

tiledlayout(1,3)

cmdurs = plasma(3) ; 

durcats = [ 0 1.4 ; 1.4 2.88 ; 2.88 10 ; 10 100 ] ; 
durcats2 = [ 0 1.4 ; 1.4 2.88 ; 2.88 5.04 ; 10 100 ] ; 

nexttile()

h = imagesc(durrmat) ; 
h.AlphaData = ~isnan(durrmat) ; 
cb = colorbar() ; 
cb.Label.String = 'edge above thr. duration' ;

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')

axis square
colormap(acton(100))

clim([0 max(durrmat,[],'all')])

nexttile()

h = imagesc(durrmat) ; 
h.AlphaData = ~isnan(durrmat) ; 
cb = colorbar() ; 
cb.Label.String = 'edge above thr. duration' ;

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')

axis square
colormap(acton(100))

clim([0 max(durrmat,[],'all')])

for idx = 1:3
    % high entries
    squrents = ~isnan(durrmat) & (durrmat >= durcats(idx,1) ) & (durrmat < durcats(idx,2)) ; 
    [u,v] = find(squrents) ; 
    n = size(durrmat,1) ; 
    for edx = 1:length(u)
    
        % imsc_addsquare(ca1,ca2,mirrorit,linecolor,linewidth)
    
        imsc_addsquare(1:n==u(edx), ...
            1:n==v(edx), ...
            0,cmdurs(idx,:),2)
    
    end
end

nt = nexttile()

h = imagesc(ones(size(durrmat))) ; 
%h.AlphaData = ~isnan(zeros(size(durrmat))) ; 
% h.AlphaData = 
cb.Label.String = 'edge above thr. duration' ;
axis square

for idx = 1:3
    % high entries
    squrents = ~isnan(durrmat) & (durrmat >= durcats2(idx,1) ) & (durrmat < durcats2(idx,2)) ; 
    [u,v] = find(squrents) ; 
    n = size(durrmat,1) ; 
    for edx = 1:length(u)
    
        % imsc_addsquare(ca1,ca2,mirrorit,linecolor,linewidth)
    
        imsc_addsquare(1:n==u(edx), ...
            1:n==v(edx), ...
            0,cmdurs(idx,:),2)
    
    end
end

colormap(nt,[0.8 0.8 0.8])

yticks(1:2:length(stimLenSweep))
yticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
ylabel('activity duration node A')

xticks(1:2:length(stimLenSweep))
xticklabels(num2str(stimLenSweep(1:2:length(stimLenSweep))'))
xlabel('activity duration node B')


%%

set(gcf,'Color','w')
orient(gcf,'landscape')
set(gcf,'Position',[0 100 800 400])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/durr_mat_w_colors.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%% now vary the length of both durations, MOREEEE

% neuralHz = 10 ;
% hrf = getcanonicalhrf_wmag(1/neuralHz,1/neuralHz) ; 
% 
% 
% fmriHz = 1/TR ; 
% 
% %spkThrSweep = ( 2:0.05:3 ) ; 
% stimLenSweep = ( 1:1:25 ) ; 
% 
% nrep=25 ; 
% res2 = zeros(length(stimLenSweep),length(stimLenSweep),nrep) ; 
% res2null = zeros(length(stimLenSweep),length(stimLenSweep),nrep) ; 
% 
% rng(42)
% for ndx = 1:nrep
% 
% disp(ndx)
% for idx = 1:length(stimLenSweep)
%     for jdx = 1:length(stimLenSweep)
% 
%         %% neural ts get noise before convolution
% 
%         % init ts of 60 sec
%         nts1 = zeros(neuralHz*60,1) ; 
%         nts2 = zeros(neuralHz*60,1) ; 
% 
%         %sig1Len = 1 ; 
% 
%         % insert the 'neural activity' at 1 sec 
%         insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
%         insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%         nts1(insertInds1) = 1 ; 
%         nts2(insertInds2) = 1 ; 
% 
%         % convolve with the hrf at neuralHz
%         [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%         [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%         % resample to fmriHZ
%         fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%         fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%         % adding noise
%         % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%         % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%         % add 2% signal noise 
%         fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*0.02)) ; 
%         fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*0.02)) ; 
% 
%         ee = prod([fts1(:) fts2(:)],2) ; 
% 
%         mm = meas_abv_thr_interp(ee,2.25).*TR ; 
%         if mm<0
%             error('neg length problem')
%         end
% 
%         res2(idx,jdx,ndx) =  mm(1) ;
% 
%         %% now the null
% 
%         nfts1 = generate_phase_surrogates(fts1) ; 
%         nfts2 = generate_phase_surrogates(fts2) ; 
% 
%         nee = prod([nfts1(:) nfts2(:)],2) ; 
% 
%         mm = meas_abv_thr_interp(nee,2.25).*TR ; 
%         mm = mm(1) ; 
%         if mm<0
%             error('neg length problem')
%         end
% 
%         res2null(idx,jdx,ndx) = mm(1) ; 
% 
%     end
% end
% 
% end
% 
% mat = mean(res2,3,'omitnan') ; 
% nmat = mean(res2null,3,'omitnan') ; 
% 
% %%
% 
% exemplarindz = [23 6; 15 16; 11 16 ] ; 
% exemplarindxcol = [0 0 1 ; 1 0 0 ; 0 1 0 ] ; 
% 
% TL = tiledlayout(2,1+size(exemplarindxcol,1),'TileIndexing','rowmajor') ; 
% 
% nexttile(TL)
% 
% %res2(isnan(res2)) = 0 ; 
% h = imagesc(mat) ; 
% h.AlphaData = ~isnan(mat) ; 
% cb = colorbar() ; 
% cb.Label.String = 'edge above thr. duration' ;
% 
% yticks(1:length(stimLenSweep))
% yticklabels(num2str(stimLenSweep'))
% ylabel('neural activity duration node A')
% 
% xticks(1:length(stimLenSweep))
% xticklabels(num2str(stimLenSweep'))
% xlabel('neural activity duration node B')
% 
% clim([0 max(mat,[],'all')])
% 
% for edx = 1:size(exemplarindxcol,1)
%     imsc_addsquare(1:length(stimLenSweep)==exemplarindz(edx,1), ...
%         1:length(stimLenSweep)== exemplarindz(edx,2), ...
%         0,exemplarindxcol(edx,:))
% end
% 
% for edx = 1:size(exemplarindxcol,1)
% 
% % get an exemplar
% 
% nreps = 50 ; 
% exres = cell(nreps,1) ; 
% 
% idx = exemplarindz(edx,1) ; 
% jdx = exemplarindz(edx,2) ; 
% 
% for ndx = 1:nreps
%     % init ts of 60 sec
%     nts1 = zeros(neuralHz*60,1) ; 
%     nts2 = zeros(neuralHz*60,1) ; 
% 
%     % insert the 'neural activity' at 1 sec 
%     sig1Len = 1 ; 
% 
%     insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
%     insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%     nts1(insertInds1) = 1 ; 
%     nts2(insertInds2) = 1 ; 
% 
%     % convolve with the hrf at neuralHz
%     [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%     [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%     % resample to fmriHZ
%     fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%     fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%     % adding noise
%     % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%     % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%     % add 2% signal noise 
%     fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*0.02)) ; 
%     fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*0.02)) ; 
% 
%     ee = prod([fts1(:) fts2(:)],2) ; 
% 
%     exres{ndx} = [fts1 fts2 ee] ; 
% 
% end
% 
% TL2 = tiledlayout(TL,3,1) ; 
% TL2.Layout.Tile = edx+1 ;
% 
% nexttile(TL2,[2 1])
% 
% timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
% wt = 1:find(timevec>45,1,'first') ; 
% 
% for ndx = 1:nreps
% 
%     plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
%     hold on
%     plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])
% 
% 
%     plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])
% 
% 
% end
% 
% text(0.5,0.4,{ ...
%     ['duration A: '  num2str(stimLenSweep(idx)) ' sec' ] ...
%     ['duration B: ' num2str(stimLenSweep(jdx)) ' sec'] ...
%     ['abv. thr. dur.: ' num2str(mat(idx,jdx)) ' sec'] }, ...
%     'units','normalized')
% 
% 
% xlim([0 45])
% ylim([-10 50])
% 
% title('example time series')
% legend({'node time series' '' 'edge time series'})
% 
% %
% 
% rng(42)
% 
% %nex = 1 ; 
% exnts = zeros(neuralHz*60,2) ; 
% 
% % insert the 'neural activity' at 1 sec 
% insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
% insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
% exnts(insertInds1,1) = 1 ; 
% exnts(insertInds2,2) = 1 ; 
% 
% 
% nexttile(TL2,[1 1])
% 
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+fliplr((0.5:2)),'Color','b')
% % plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
% ylim([0 2.5])
% xlim([0 45])
% 
% yticks([.5 1.5])
% yticklabels({'B' 'A'})
% 
% title('example boxcar')
% xlabel('time (sec)')
% 
% end
% 
% %
% 
% nexttile(TL)
% 
% %res2(isnan(res2)) = 0 ; 
% 
% h = imagesc(nmat) ; 
% h.AlphaData = ~isnan(nmat) ; 
% cb = colorbar() ; 
% cb.Label.String = 'edge above thr. duration' ;
% 
% yticks(1:length(stimLenSweep))
% yticklabels(num2str(stimLenSweep'))
% ylabel('neural activity duration node A')
% 
% xticks(1:length(stimLenSweep))
% xticklabels(num2str(stimLenSweep'))
% xlabel('neural activity duration node B')
% 
% clim([0 max(mat,[],'all')])
% 
% for edx = 1:size(exemplarindxcol,1)
%     imsc_addsquare(1:length(stimLenSweep)==exemplarindz(edx,1), ...
%         1:length(stimLenSweep)== exemplarindz(edx,2), ...
%         0,exemplarindxcol(edx,:))
% end
% 
% for edx = 1:size(exemplarindxcol,1)
% 
% % get an exemplar
% 
% nreps = 50 ; 
% exres = cell(nreps,1) ; 
% 
% idx = exemplarindz(edx,1) ; 
% jdx = exemplarindz(edx,2) ; 
% 
% for ndx = 1:nreps
%     % init ts of 60 sec
%     nts1 = zeros(neuralHz*60,1) ; 
%     nts2 = zeros(neuralHz*60,1) ; 
% 
%     % insert the 'neural activity' at 1 sec 
%     sig1Len = 1 ; 
% 
%     insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
%     insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%     nts1(insertInds1) = 1 ; 
%     nts2(insertInds2) = 1 ; 
% 
%     % convolve with the hrf at neuralHz
%     [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%     [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%     % resample to fmriHZ
%     fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%     fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%     % adding noise
%     % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%     % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%     % add 2% signal noise 
%     fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*0.02)) ; 
%     fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*0.02)) ; 
% 
%     nfts1 = generate_phase_surrogates(fts1);
%     nfts2 = generate_phase_surrogates(fts2);
% 
%     nee = prod([nfts1(:) nfts2(:)],2) ; 
% 
%     exres{ndx} = [nfts1 nfts2 nee] ; 
% 
% end
% 
% TL2 = tiledlayout(TL,3,1) ; 
% TL2.Layout.Tile = edx+5 ;
% 
% nexttile(TL2,[2 1])
% 
% timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
% wt = 1:find(timevec>45,1,'first') ; 
% 
% for ndx = 1:nreps
% 
%     plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
%     hold on
%     plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])
% 
% 
%     plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])
% 
% 
% end
% 
% hr = (nrep-sum(isnan(res2null),3))./100 ; 
% 
% text(0.1,0.7,{ ...
%     ['duration A: '  num2str(stimLenSweep(idx)) ' sec' ] ...
%     ['duration B: ' num2str(stimLenSweep(jdx)) ' sec'] ...
%     ['abv. thr. dur.: ' num2str(nmat(idx,jdx)) ' sec'] ...
%     ['hit rate: ' num2str(hr(idx,jdx)) '%' ] }, ...
%     'units','normalized')
% 
% 
% xlim([0 45])
% ylim([-2 10])
% 
% title('example time series')
% legend({'node time series' '' 'edge time series'})
% 
% %
% 
% rng(42)
% 
% %nex = 1 ; 
% exnts = zeros(neuralHz*60,2) ; 
% 
% % insert the 'neural activity' at 1 sec 
% insertInds1 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(idx)*neuralHz)-1) ;
% insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
% exnts(insertInds1,1) = 1 ; 
% exnts(insertInds2,2) = 1 ; 
% 
% 
% nexttile(TL2,[1 1])
% 
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+fliplr((0.5:2)),'Color','b')
% % plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
% ylim([0 2.5])
% xlim([0 45])
% 
% yticks([.5 1.5])
% yticklabels({'B' 'A'})
% 
% title('example boxcar')
% xlabel('time (sec)')
% 
% end
% 
% %%
% 
% set(gcf,'Color','w')
% orient(gcf,'landscape')
% set(gcf,'Position',[0 100 1400 700])
% 
% out_figdir = [ './reports/figures/supp/' ]
% mkdir(out_figdir)
% filename = [out_figdir '/gridBIGGER_null_mats.pdf' ] ; 
% print(filename,'-dpdf','-vector')
% close(gcf)

%% 

% neuralHz = 10 ; 
% 
% poissonSpkGen = @(firingrateHz,lenSec) rand(1,floor(lenSec*neuralHz)) < (firingrateHz*1/neuralHz) ; 
% regularSpkGen = @(firingrateHz,lenSec) ismember(1:lenSec*neuralHz,floor(linspace(1,lenSec*neuralHz,floor(firingrateHz*lenSec)))) ; 
% irregularSpkGen = @(firingrateHz,lenSec,poissShiftVar) ...
%     ismember(1:lenSec*neuralHz, arrayfun(@(x_) x_ + (poissrnd(poissShiftVar)*randsample([-1,1],1)) , find(regularSpkGen(firingrateHz,lenSec)) )) ; 
% boxCarSpkGen = @(firingrateHz,lenSec) ones(floor(lenSec*neuralHz),1) ; 
% 
% hrf = getcanonicalhrf_wmag(1/neuralHz,1/neuralHz) ; 
% 
% useGenFunc = boxCarSpkGen ; 
% 
% %
% 
% fmriHz = 1/TR ; 
% 
% spkHzSweep = ( 1:0.5:10 ) ; 
% stimLenSweep = ( 1:0.5:10 ) ; 
% 
% nrep=50 ; 
% res1 = zeros(length(spkHzSweep),length(stimLenSweep),nrep) ; 
% res2 = zeros(length(spkHzSweep),length(stimLenSweep),nrep) ; 
% 
% rng(42)
% for ndx = 1:nrep
% 
% disp(ndx)
% for idx = 1:length(spkHzSweep)
%     for jdx = 1:length(stimLenSweep)
% 
%         %% neural ts get noise before convolution
% 
%         % init ts of 60 sec
%         nts1 = zeros(neuralHz*60,1) ; 
%         nts2 = zeros(neuralHz*60,1) ; 
% 
%         sig1Len = 1 ; 
% 
%         % insert the 'neural activity' at 1 sec 
%         insertInds1 = (2*neuralHz):((2*neuralHz)+(sig1Len*neuralHz)-1) ;
%         insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%         nts1(insertInds1) = useGenFunc(spkHzSweep(idx),sig1Len) ; 
%         nts2(insertInds2) = useGenFunc(spkHzSweep(idx),stimLenSweep(jdx)) ; 
% 
%         % convolve with the hrf at neuralHz
%         [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%         [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%         % resample to fmriHZ
%         fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%         fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%         % adding noise
%         % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%         % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%         % add 2% signal noise 
%         fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*0.02)) ; 
%         fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*0.02)) ; 
% 
%         ee = prod([fts1(:) fts2(:)],2) ; 
% 
%         mm = meas_abv_thr_interp(ee,2.25).*TR ; 
%         if mm<0
%             error('neg length problem')
%         end
% 
%         res2(idx,jdx,ndx) =  mm(1) ;
% 
%     end
% end
% 
% end
% 
% mat = mean(res2,3,'omitnan') ; 
% 
% %%
% 
% % set noise to %2 of signal
% % do boxcar
% % take a few thresholds 
% 
% exemplarindz = [13 9; 12 16; 4 19 ] ; 
% exemplarindxcol = [0 0 1 ; 1 0 0 ; 0 1 0 ] ; 
% 
% TL = tiledlayout(1,1+size(exemplarindxcol,1)) ; 
% 
% nexttile(TL)
% 
% 
% %res2(isnan(res2)) = 0 ; 
% h = imagesc(mat) ; 
% h.AlphaData = ~isnan(mat) ; 
% cb = colorbar() ; 
% cb.Label.String = 'edge above thr. duration' ;
% 
% yticks(1:length(spkHzSweep))
% yticklabels(num2str(spkHzSweep'))
% ylabel('neural spk. Hz')
% 
% xticks(1:length(stimLenSweep))
% xticklabels(num2str(stimLenSweep'))
% xlabel('stim duration')
% 
% 
% clim([0 max(mat,[],'all')])
% 
% for edx = 1:size(exemplarindxcol,1)
%     imsc_addsquare(1:length(spkHzSweep)==exemplarindz(edx,1), ...
%         1:length(stimLenSweep)== exemplarindz(edx,2), ...
%         0,exemplarindxcol(edx,:))
% end
% 
% for edx = 1:size(exemplarindxcol,1)
% 
% % get an exemplar
% 
% nreps = 50 ; 
% exres = cell(nreps,1) ; 
% 
% idx = exemplarindz(edx,1) ; 
% jdx = exemplarindz(edx,2) ; 
% 
% for ndx = 1:nreps
%     % init ts of 60 sec
%     nts1 = zeros(neuralHz*60,1) ; 
%     nts2 = zeros(neuralHz*60,1) ; 
% 
%     % insert the 'neural activity' at 1 sec 
%     sig1Len = 1 ; 
% 
%     insertInds1 = (2*neuralHz):((2*neuralHz)+(sig1Len*neuralHz)-1) ;
%     insertInds2 = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
%     nts1(insertInds1) = useGenFunc(spkHzSweep(idx),sig1Len) ; 
%     nts2(insertInds2) = useGenFunc(spkHzSweep(idx),stimLenSweep(jdx)) ; 
% 
%     % convolve with the hrf at neuralHz
%     [~,cnts1] = pad_conv(nts1,hrf,100) ; 
%     [~,cnts2] = pad_conv(nts2,hrf,100) ; 
% 
%     % resample to fmriHZ
%     fts1 = resampsig1d(cnts1,neuralHz,fmriHz)  ; 
%     fts2 = resampsig1d(cnts2,neuralHz,fmriHz)  ; 
% 
%     % adding noise
%     % fts1 = fts1 + normrnd(0,0.05,size(fts1)) ; 
%     % fts2 = fts2 + normrnd(0,0.05,size(fts2)) ;
% 
%     % add 2% signal noise 
%     fts1 = fts1 + ((rand(size(fts1))-0.5).*(max(fts1).*0.02)) ; 
%     fts2 = fts2 + ((rand(size(fts2))-0.5).*(max(fts2).*0.02)) ; 
% 
%     ee = prod([fts1(:) fts2(:)],2) ; 
% 
%     exres{ndx} = [fts1 fts2 ee] ; 
% 
% end
% 
% %
% 
% TL2 = tiledlayout(TL,3,1) ; 
% TL2.Layout.Tile = edx+1 ;
% 
% nexttile(TL2,[2 1])
% 
% timevec = (0:ceil(60*fmriHz)-1).*(1/fmriHz) ; 
% wt = 1:find(timevec>45,1,'first') ; 
% 
% for ndx = 1:nreps
% 
%     plot(timevec(wt)',exres{ndx}(wt,1),'Color',[0.5 0.5 0.5 0.1])
%     hold on
%     plot(timevec(wt)',exres{ndx}(wt,2),'Color',[0.5 0.5 0.5 0.1])
% 
% 
%     plot(timevec(wt)',exres{ndx}(wt,3),'Color',[exemplarindxcol(edx,:) 0.3])
% 
% 
% end
% 
% text(0.6,0.6,['spike Hz: ' num2str(spkHzSweep(idx))],'units','normalized')
% text(0.6,0.5,['spiking duration: ' num2str(stimLenSweep(jdx)) ' sec.'],'units','normalized')
% text(0.6,0.4,['abv thr len: ' num2str(mat(idx,jdx)) ' sec.'],'units','normalized')
% 
% xlim([0 45])
% ylim([-2 20])
% 
% title('example time series')
% legend({'node time series' '' 'edge time series'})
% 
% %
% 
% rng(42)
% 
% nex = 1 ; 
% exnts = zeros(neuralHz*60,nex) ; 
% 
% for tdx = 1:nex
% 
% % insert the 'neural activity' at 1 sec 
% insertInds = (2*neuralHz):((2*neuralHz)+(stimLenSweep(jdx)*neuralHz)-1) ;
% 
% exnts(insertInds,tdx) = useGenFunc(spkHzSweep(idx),stimLenSweep(jdx)) ; 
% 
% end
% 
% nexttile(TL2,[1 1])
% 
% plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts+0.5),'Color','b')
% % plot((0:1:(neuralHz*60)-1)/neuralHz,(exnts.*0.5)+(0.5:5),'Color','b')
% ylim([0 5.5])
% xlim([0 45])
% 
% title('example neural spikes')
% 
% end