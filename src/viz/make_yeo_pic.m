neword = [ 12 13 11 16 17 15 5 6 10 9 7 8 3 4 14 1 2 ];

ccc = [ 

    0.4688    0.0703    0.5234
    0.9961         0         0
    0.2734    0.5078    0.7031
    0.1641    0.7969    0.6406
    0.2891    0.6055    0.2344
         0    0.4609    0.0547
    0.7656    0.2266    0.9766
    0.9961    0.5938    0.8320
    0.8594    0.9688    0.6406
    0.4766    0.5273    0.1953
    0.4648    0.5469    0.6875
    0.8984    0.5781    0.1328
    0.5273    0.1953    0.2891
    0.0469    0.1875    0.9961
         0         0    0.5078
    0.9961    0.9961         0
    0.8008    0.2422    0.3047
 ] ; 

ca = [

    16
    16
    16
    16
    16
    16
    17
    17
    17
    17
    17
    17
    13
    13
    13
    13
    13
    13
    13
    13
    14
    14
    14
    14
    14
    14
    14
    14
     7
     7
     7
     7
     7
     7
     8
     8
     8
     8
     8
    11
    11
    11
    11
    11
    11
    11
    12
    12
    12
    12
     9
     9
    10
    10
    10
    10
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     2
     2
     2
     2
     2
     3
     3
     3
     4
     4
     4
     4
     4
     4
     4
     4
     5
     5
     5
     5
     5
     5
     5
     5
     5
     5
     5
     5
     5
     6
     6
     6
    15
    15
    16
    16
    16
    16
    16
    16
    17
    17
    17
    17
    17
    17
    13
    13
    13
    13
    13
    13
    13
    13
    13
    13
    13
    14
    14
    14
    14
    14
    14
    14
     7
     7
     7
     7
     7
     7
     8
     8
     8
     8
     8
    11
    11
    11
    11
    11
    11
    11
    11
    11
    12
    12
    12
    12
    12
    12
     9
     9
     9
     9
    10
    10
    10
    10
     1
     1
     1
     1
     1
     1
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     3
     3
     3
     4
     4
     4
     4
     4
     4
     5
     5
     5
     5
     6
     6
     6
    15
    15
    15
    15 ] ; 

parc_plot(surfss,annotm,'schaefer200-yeo17', ca,...
    'cmap',ccc(neword,:),...
    'viewcMap',0,'newFig',0,'viewStr','all')

set(gcf,'Position',[100 100 800 600])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/yeo_comms.pdf' ] ; 
print(filename,'-dpdf')
close(gcf)

%% and the colorbar

cc = ccc(neword,:) ; 
pp = parc.names(neword) ; 
[~,si] = sort(pp) ; 

imagesc((1:17)')
colormap(cc)
xticks([])
yticks(1:1:17)
yticklabels(pp(si))


yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'none'; 

set(gca,'TickLength',[0 0])
set(gcf,'Position',[100 100 125 400])

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/yeo_comms_names.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)