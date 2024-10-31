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

tiledlayout(1,3) ;

hrfamps = [1.7 1.85 2 ] ; 

for tdx = 1:3

    hrfamp = hrfamps(tdx) ; 
    hrf = getcanonicalhrf(.1,1)*hrfamp ; 
    
    hrf2_onsets = 0:0.5:5 ; 
    hrf2_lengths = 0.1:0.5:5 ; 
    % magDelta = 1:0.1:4 ; 
    
    res = zeros(length(hrf2_onsets),length(hrf2_lengths)) ; 
    
    for idx = 1:length(hrf2_onsets)
        for jdx = 1:length(hrf2_lengths)
    
            hrf1 = hrf.*1 ; 
            hrf2 = getcanonicalhrf(hrf2_lengths(jdx),1)*hrfamp ; 
    
            h1 = [ hrf1 ] ;
            
            h2 = interp1((0:length(hrf2)-1)+hrf2_onsets(idx),hrf2,(0:length(hrf2)-1),"linear") ; 
            h2 = h2(1:length(h1)) ; 
            h2(isnan(h2)) = 0 ; 
    
            ee = prod([h1(:) h2(:)],2) ; 
    
            res(idx,jdx) = meas_abv_thr_interp(ee,2.25) ; 
    
        end
    end

    nexttile()

    imagesc(res)

end

%% 



tr = 100 ; 

poissonSpkGen = @(firingrateHz,lenSec) rand(1,floor(lenSec*tr)) < (firingrateHz*1/tr) ; 
regularSpkGen = @(firingrateHz,lenSec) ismember(1:lenSec*tr,floor(linspace(1,lenSec*tr,floor(firingrateHz*lenSec)))) ; 
irregularSpkGen = @(firingrateHz,lenSec,poissVar) ...
    ismember(1:lenSec*tr, arrayfun(@(x_) x_ + (poissrnd(poissVar)*randsample([-1,1],1)) , find(regularSpkGen(firingrateHz,lenSec)) )) ; 


hrf = getcanonicalhrf(.1,1/tr) ; 

neuralhz = [  ]



 ismember(1:lenSec*tr,floor(linspace(1,lenSec*tr,floor(firingrateHz*lenSec))))