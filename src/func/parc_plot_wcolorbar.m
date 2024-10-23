function [pp] = parc_plot_wcolorbar(dat,surfss,annotm,valrange,cm,figsz)

if nargin < 3
    error('need at least 3 args')
end

if nargin < 4
    valrange = [min(dat) max(dat)] ; 
end

if nargin < 5
    cm = parula() ; 
end

if nargin < 6
    figsz = [100 100 600 1000] ; 
end

TL = tiledlayout(2,1) ;
nt1 = nexttile() ;

pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', dat ,...
    'valRange',valrange,...
    'cmap',cm, ...
    'viewcMap',0,'newFig',0,'viewStr','all',...
    'parenth',TL) ;
pp.Layout = nt1.Layout ; 

nt2 = nexttile(TL) ; 

hh = imagesc(dat) ;
cb = colorbar() ;
clim(valrange)
hh.Visible = 'off' ;
hh.Parent.Visible = 'off' ; 
cb.Location = "north" ; 
cl = clim() ; 
cb.Ticks = linspace(cl(1),cl(2),5) ; 
cb.TickLabels = strtrim(cellstr(num2str(cellfun(@(x_) round(str2num(x_),2) , cb.TickLabels)))) ; 
cm = colormap() ; 
colormap(nt2,cm(2:end,:))

TL.TileSpacing = 'tight' ; 

set(gcf,'Position',figsz)
set(gcf,'Color','w')

