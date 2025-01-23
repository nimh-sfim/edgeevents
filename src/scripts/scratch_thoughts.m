
a = tv(meanc) ; 
b = tv(meanlen) ; 
c = tv(ietmean) ; 


AA = discretize(a,prctile(a,0:100))./100 ; 
BB = discretize(b,prctile(b,0:100))./100 ; 
CC = discretize(c,prctile(c,0:100))./100 ; 

scatter(b,c,30,a,'filled')

%%

ca = parc.ca(1:200) ; 
res = zeros((finfo.nnodes*(finfo.nnodes-1))/2,17) ; 

for idx = 1:17
    tmp = (ca==idx).*(ca==idx)' ; 
    res(:,idx) = tv(tmp) ; 

end

%%

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

newcol = ccc(neword,:) ; 

scatter(a,b,30,'filled','MarkerFaceColor',[0.8 0.8 0.8]) 

for idx = 1:17

%ii = sum(res,2)==0 ; 

hold on
iii = ~~res(:,idx) ; 
ai = a(iii) ; bi = b(iii) ; 
scatter(a(iii),b(iii),40,'filled','MarkerFaceColor',newcol(idx,:))

k = convhull(ai,bi) ;
fill(ai(k),bi(k),newcol(idx,:),'FaceAlpha',0.3)

% title(idx)
% waitforbuttonpress
% clf
end

%%

xdat = tv(meanc) ; 
ydat = tv(meanlen) ; 

xbins = prctile(xdat,0:1:100) ; % just use simple percentiles, all bins have
                               % same numeber of edges
[eVar,xbinmid,binmean] = norm_bin_model(xdat,xbins,ydat) ; 