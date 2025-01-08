%% make some basic viz stuff

idx = 42 ; 
ts = datStr(idx).ts(:,1:200) ; 
ets = get_ets(ts) ; 
fc = corr(ts) ; 

[u,v] = find(triu(ones(200),1)) ; 

%%

timevec = 1:150 ; 
trvec = timevec.*0.72 ; 

%%

t1 = zscore(ts(:,u(15))) ; 
t2 = zscore(ts(:,v(15))) ; 

ets1 = prod([ zscore(t1) zscore(t2)],2) ; 

[c1,c2] = count_spks(ets1,2.25) ;
[s1,s2] = spk_lenmat(ets1>2.25) ;

cc2 = nan(size(c2)) ; 
cc2(c2) = 1 ; 

cc3 = nan(size(cc2)) ; 
cc3(ets1>2.25) = 1 ; 

%%


tiledlayout(5,1)
nexttile()

plot(trvec,t1(timevec),'.-')
ylim([-5 5])
yline(0)
xlim([min(trvec) max(trvec)])

nexttile()
plot(trvec,t2(timevec),'.-')
ylim([-5 5])
yline(0)
xlim([min(trvec) max(trvec)])

nexttile()
plot(trvec,ets1(timevec),'.-')
ylim([-5 10])
yline(0)
xlim([min(trvec) max(trvec)])

nexttile()
axis off
nexttile()
plot(trvec,ets1(timevec),'.-')
ylim([-5 10])
yline(0)

hold on
% plot(trvec',cc2(timevec).*9.5,'*')

numlabtp = find(~isnan(cc2(timevec))) ; 
for idx = numlabtp'
    text(idx.*0.72,15,...
        [ num2str(s1(idx).*0.72) ],...
        'Color','r','HorizontalAlignment','center','VerticalAlignment','cap','FontSize',11,...
        'Rotation',45)
    line([idx.*0.72 idx.*0.72],[2.25 12],'Color','r')

end

plot(trvec',cc3(timevec).*ets1(timevec),'.','MarkerSize',10,'Color','r')

yline(2.25,'Color','r','LineStyle','--')
xlim([min(trvec) max(trvec)])

set(gcf,'Position',[100 100 800 600])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figInfo/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/example_ts_w_spike.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)

%%


minifc = fc(1:20,1:20) ; 

imsc_grid_comm(minifc,1:20,0.5,[0.5 0.5 0.5],[0.2 .2 .2])
colormap(parula)
clim([-0.5 1])

axis off
axis square

out_figdir = [ './reports/figures/figInfo/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/example_minimat.pdf' ] ; 
print(filename,'-dpdf','-bestfit')
close(gcf)


%%

[vec,~] = eigs(double(meanfc),1) ; 

examplenodes = normalize(vec,'range',[0 0.3]) ; 
examplenodes(((1:200)==151)|((1:200)==155)) = 1 ; 

% examplenodes = zeros(200,1) ; 
% examplenodes = 1:200 ; 



%151
%155

% ((1:200)==151)|((1:200)==155)

p = parc_plot(surfss,annotm,'schaefer200-yeo17',...
    examplenodes,...
            'cmap', [parula(100)] ,...
        'viewcMap',0,'newFig',0,'viewStr','rh:lat') ;

        % 'cmap',[0.5 0.5 0.5 ; [102 170 215]./255 ] ,...


set(gcf,'Position',[100 100 800 600])
set(gcf,'Color','w')

out_figdir = [ './reports/figures/figInfo/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/example_surface.png' ] ; 
print(filename,'-dpng')
close(gcf)

%%


idx = 42 ; 
ts = datStr(idx).ts(:,1:200) ; 
ets = get_ets(ts) ; 
fc = corr(ts) ; 

[u,v] = find(triu(ones(200),1)) ; 

timevec = 1:150 ; 
trvec = timevec.*0.72 ; 

%%
 
miniets = flipud(ets(timevec,442:452))+((1:11)) ; 

% dat = miniets ; 
hold off
for idx = 1:size(miniets,2)
    % plot(-42:1:42,dat(:,idx),'LineWidth',2)
    x = 1:size(miniets,1) ; 
    y = miniets(:,idx)' ; 
    z = zeros(size(y)) ; 
    % c = (miniets(:,idx))' ; 
    c = single(miniets(:,idx)>(2.25+idx))' ; 
    surface([x;x],[y;y],[z;z],[c;c],...
        'FaceColor','no','EdgeColor','interp','LineWidth',2,'EdgeAlpha',1)
end
hold off
colormap([0.5 0.5 0.5; 1 0 0 ])

% for idx = 1:11
%     plot(miniets(:,idx),'Color',[0.2 0.2 0.2])
%     hold on
% 
%     ab = miniets(:,idx)>(2.25+idx) ; 
%     tvec = nan(150,1) ; 
%     tvec(ab) = miniets(ab,idx) ; 
%     plot(timevec,tvec,'Color','r','LineWidth',2,'LineStyle',':')
% 
% end

xticks([])
yticks([])

out_figdir = [ './reports/figures/figInfo/' ]
orient(gcf,'landscape')
mkdir(out_figdir)
filename = [out_figdir '/spk_later.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% 80 percentile

m1 = prctile(tv(spike_conn.subset1.long),79.9) ; 
m2 = prctile(tv(spike_conn.subset1.long),80.1) ; 

m3 = (spike_conn.subset1.long>=m1) & (spike_conn.subset1.long<=m2)
