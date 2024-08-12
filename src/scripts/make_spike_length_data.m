%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

%%

SPK_THR = 2 ; 

for idx = 1:NSUBS

    disp(idx)

    filename = [DD.PROC '/' datStr(idx).sub '_' OUTSTR '_spike_len.mat'] ; 

    if isfile(filename)
        disp(['already finsiehd sub: ' datStr(idx).sub ])
        continue
    end

    tmpts = datStr(idx).ts ; 
    tmpets = get_ets(tmpts) ; 

    abv_thr = tmpets > SPK_THR ; 
    [spike_len_mat,spike_len_cell] = spk_lenmat(abv_thr) ; 

    save([DD.PROC '/' datStr(idx).sub '_' OUTSTR '_spike_len.mat'],...
        'spike_len_cell','spike_len_mat','-v7.3')

end

%% now read in spike lengths to make a histogram of lengths

subsets = {'subset1' 'subset2'} ; 

for sdx = subsets

    spike_lengths.(sdx{1}) = cell(length(sublist.(sdx{1})),1) ; 
    
    for idx = 1:length(sublist.(sdx{1}))
    
        disp(idx)
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' datStr(sind).sub '_' OUTSTR '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_cell') ; 
    
        spike_lengths.(sdx{1}){idx} = readdat.spike_len_cell  ; 
    
    end

end

%%

for sdx = subsets
    tmp = cell2mat(arrayfun(@(i_) int32(cell2mat(spike_lengths.(sdx{1}){i_}')),1:length(spike_lengths.subset1)','UniformOutput',false)') ;
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

%% 

% look at the subjects 
imagesc(spike_prct.subset1.prct(:,fliplr(sortedInd(mean(spike_prct.subset1.prct(95:100,:))))))

%% find way to naturally bin it

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

%%

e1 = [ 1 2 4 6 maxspk ] ; 
e2 = [ 1 2 3 4 5 maxspk ] ; 



% 
% %%
% 
% spike_prct.subset1.agg.by5 = prctile(all_spike_lengths, 5:5:100 ) ; 
% spike_prct.subset1.agg.bylog10 = prctile(all_spike_lengths, logspace(log10(1),log10(100),20) ) ; 
% spike_prct.subset1.agg.byexp2 = prctile(all_spike_lengths, linspace(1,10,20).^2 ) ; 
% 
% %% look at spike length for each pers
% 
% spike_prct.subset1.by5 = nan(20,length(sublist.subset1)) ; 
% spike_prct.subset1.byexp2 = nan(20,length(sublist.subset1)) ; 
% 
% for idx = 1:length(sublist.subset1)
% 
%     disp(idx)
% 
%     tmp = int32(cell2mat(spike_lengths.subset1{idx}')) ; 
% 
%     spike_prct.subset1.by5(:,idx) = prctile(tmp, 5:5:100 ) ; 
%     spike_prct.subset1.byexp2(:,idx) = prctile(tmp, linspace(1,10,20).^2 ) ; 
% end
% 
% %% divide up 
% 
% maxspk = max(spike_prct.subset1.agg.by5) ; 
% thrs = [ 2 3 4 maxspk ] ; 
% prct_inc = 5:5:100 ; 
% 
% prct_vals = nan(length(thrs),1) ; 
% for idx = 1:length(thrs)
%     prct_vals(idx) = prct_inc(find(spike_prct.subset1.agg.by5==thrs(idx),1,'last')) ; 
% end
% 
% %%

%% 

% tiledlayout

dd = discretize(spike_lengths.all.subset1,e1) ; 
ca = parula(max(dd)) ;

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:maxspk,FaceColor=ca(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    sum(dd==idx)./length(dd)
    hold on

end
hold off

%%

dd = discretize(spike_lengths.all.subset1,[1 5 9 maxspk]) ; 
ca = parula(max(dd)) ;

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:maxspk,FaceColor=ca(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    sum(dd==idx)./length(dd)
    hold on

end
hold off



%% 

dd = discretize(spike_lengths.all.subset1,e2) ; 
ca = parula(max(dd)) ;

for idx = 1:max(dd)

    histogram(spike_lengths.all.subset1(dd==idx), 0:1:maxspk,FaceColor=ca(idx,:)) ;
    % [min(spike_lengths.all.subset1(dd==idx)) max(spike_lengths.all.subset1(dd==idx))]
    sum(dd==idx)./length(dd)
    hold on

end
hold off

