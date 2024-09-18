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

fullunroll = @(x_) x_(:) ; 
subsets = {'subset1' 'subset2'} ; 

refparc = parc.ca(1:200) ; 

%% now read in spike lengths to make a histogram of lengths

spike_hubs = struct() ; 

for sdx = subsets

    for hdx = 1:3
        spike_hubs.(sdx{1}).(spklen_names{hdx}) = zeros(200,length(sublist.(sdx{1}))) ; 
    end

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        origfc = corr(datStr(sind).ts(:,1:200)) ; 

        m = single(percthr_trimlow_und(origfc)) ; 
        m(1:finfo.nnodes+1:end) = 0 ; 

        ss = strengths_und(m.*origfc) ; 
        orighub = ss>=prctile(ss,80) ; 

        %% lets figure out if the hubs are the same! 

        for hdx = 1:3

            mat = mksq(mean(dd==hdx)) ; 

            % [h,hh,hhh] = get_hub_score_wei_und(abs(mat),'gollo2018','zscore',.20) ; 
            % hh = participation_coef(abs(mat.*m),refparc) ; 
    
            matstr = strengths_und(mat) ; 
            
            spike_hubs.(sdx{1}).(spklen_names{hdx})(:,idx) = matstr>=prctile(matstr,80) ; 
        end

    end
end

%%

figure
cortexplot(mean(spike_hubs.subset2.(spklen_names{1}),2))

%% divide out by yeo sys

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 

remap_nodes = remaplabs(parc.ca(1:200),g1sort,1:17) ; 
[datplot_ca,datplot_sinds] = sort(remap_nodes) ;  
datplot_lab_sinds = remaplabs(1:17,g1sort,1:17) ; 

cmap = get_nice_yeo_cmap('grad1') ; 

clf
TL = tiledlayout(2,3) ; 

for idx = 1:3

    nt1 = nexttile(TL)
    
    tmp = mean(spike_hubs.subset1.(spklen_names{idx}),2) ;
    
    pp =  parc_plot(surfss,annotm,'schaefer200-yeo17', tmp ,...
        'valRange',[0 1],...
        'cmap',parula(100), ...
        'viewcMap',0,'newFig',0,'viewStr','all',...
        'parenth',TL)
    pp.Layout = nt1.Layout ; 


end

for idx = 4:6

    nexttile(TL)
    
    dat = mean(spike_hubs.subset1.(spklen_names{idx-3}),2) ; 
    
    
    fcn_boxpts(dat(datplot_sinds),...
        datplot_ca,cmap(datplot_lab_sinds,:),...
        0,parc.names(g1sort))
    
    ylim([0 1])

end

% clf
% fcn_boxpts(dat,parc.ca(1:200),get_nice_yeo_cmap('grad1'),0,parc.names(1:17))

% nothing burger 
