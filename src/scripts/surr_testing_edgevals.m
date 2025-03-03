%% start
%% TODO ADAPT THIS FOR SIG TESTING ALL EDGES 

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

g1sort = [ 13 17 14 8 16 11 7 15  12 10  1 3 6 9 2  4 5 ] ; 
remap_labs = remaplabs(parc.ca(1:200),g1sort,1:17) ;

%% load the spike conn

filename = [ DD.PROC '/spk_conn_avg_' OUTSTR '.mat' ] ; 
load(filename)

lennames = {'short' 'inter' 'long'} ; 

%% do some stat testing

% load up the nulls
filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spkcount.mat'] ; 
load(filename)

%% boots

nboot = 5000 ; 

nsubs = length(sublist.all) ; 
subsetsz = length(sublist.subset1) ; 

nperms = length(simmat.keepcov) ; 
rng(42)

surrA = struct() ;
surrB = struct() ;

% initialize
for sdx = 1:3
    surrA.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot); 
    surrB.(lennames{sdx}) = zeros(finfo.nnodes,finfo.nnodes,nboot) ; 
end

for idx = 1:nboot
    
    disp(idx)

    % pickout the surr data to average
    bootinds = randsample(1:nperms,subsetsz,true) ; 

    tmpA = zeros(3,finfo.nnodes*(finfo.nnodes-1)/2) ; 
    tmpB = zeros(3,finfo.nnodes*(finfo.nnodes-1)/2)  ; 

    for bdx = bootinds 

        aa = arrayfun(@(i_)simmat.nocov{bdx}(i_,:),1:3,'UniformOutput',false) ; 
        bb = arrayfun(@(i_)simmat.keepcov{bdx}(i_,:),1:3,'UniformOutput',false) ; 

        for sdx = 1:length(lennames)
            tmpA(sdx,:) = tmpA(sdx,:) + aa{sdx}  ; 
            tmpB(sdx,:) = tmpB(sdx,:) + bb{sdx}  ; 
        end

    end
    
    % if any(any(mksq(tmpB(1,:))<mksq(tmpB(2,:))))
    %     error('whats up here')
    % end

    % make mean
    for sdx = 1:2
        tmpA(sdx,:) = tmpA(sdx,:) ./ subsetsz  ; 
        tmpB(sdx,:) = tmpB(sdx,:) ./ subsetsz  ; 
    end

    % and record for this boot
    for sdx = 1:length(lennames)
        surrA.(lennames{sdx})(:,:,idx) = mksq(tmpA(sdx,:)) ; 
        surrB.(lennames{sdx})(:,:,idx) = mksq(tmpB(sdx,:)) ; 
    end

end 

%%

% A: no cov
% B: keep cov

sig = struct() ; 

for sdx = 1:length(lennames)

    disp(sdx)

    % more count

    Amore = sum(surrA.(lennames{sdx})>=spike_conn.subset1.(lennames{sdx}),3) ; 
    Bmore = sum(surrB.(lennames{sdx})>=spike_conn.subset1.(lennames{sdx}),3) ; 

    sig.A.more.mat.(lennames{sdx}) = Amore ; 
    sig.B.more.mat.(lennames{sdx}) = Bmore ; 

    sig.A.more.fdr.(lennames{sdx}) = mksq( fdr_bh(tv((Amore+1)./(nboot+1))) ) ; 
    sig.B.more.fdr.(lennames{sdx}) = mksq( fdr_bh(tv((Bmore+1)./(nboot+1))) ) ; 

    % less

    Aless = sum(surrA.(lennames{sdx})<=spike_conn.subset1.(lennames{sdx}),3) ; 
    Bless = sum(surrB.(lennames{sdx})<=spike_conn.subset1.(lennames{sdx}),3) ; 

    sig.A.less.mat.(lennames{sdx}) = Aless ; 
    sig.B.less.mat.(lennames{sdx}) = Bless ; 

    sig.A.less.fdr.(lennames{sdx}) = mksq( fdr_bh(tv((Aless+1)./(nboot+1))) ) ; 
    sig.B.less.fdr.(lennames{sdx}) = mksq( fdr_bh(tv((Bless+1)./(nboot+1))) ) ; 

end

%%

tiledlayout(2,4,'TileIndexing','columnmajor')

for sdx = 1:3

    nexttile()

    empmat = spike_conn.subset1.(lennames{sdx}) ; 

    if sdx == 1
        extras = {1 '' [0.1 0.1 0.1 0.5] parc.names(g1sort)} ; 
    else
        extras = {1 '' [0.1 0.1 0.1 0.5]} ; 
    end

    sigmat = sig.B.more.fdr.(lennames{sdx}) ; 
    h = imsc_grid_comm(empmat.*sigmat ,remap_labs,extras{:}) ; 
    h.AlphaData = sigmat ~= 0 ; 
    xticks('')
    axis square

    if sdx ~= 1
        yticks('') ; 
    end

    xlabel('sig. more than surrogate')
    title(lennames{sdx})

    nexttile()

    sigmat = sig.B.less.fdr.(lennames{sdx}) ; 
    h = imsc_grid_comm(empmat.*sigmat ,remap_labs,extras{:}) ; 
    h.AlphaData = sigmat ~= 0 ; 
    xticks('')
    axis square

    xlabel('sig. less than surrogate')
    
    if sdx ~= 1
        yticks('') ; 
    end

end

set(gcf,'Position',[100 100 800 500])
set(gcf,'Color','w')
orient(gcf,'landscape')

% add the ultra-long!
filename = [ DD.PROC '/spk_longestedge_sigmask_' OUTSTR '.mat' ] ; 
ll = load(filename) ; 

filename = [ DD.PROC '/spk_conn_ex_avg_' OUTSTR '.mat' ] ; 
mm = load(filename) ; 

nexttile()

empmat = mm.spike_conn_ex.subset1.longest ; 

if sdx == 1
    extras = {1 '' [0.1 0.1 0.1 0.5] parc.names(g1sort)} ; 
else
    extras = {1 '' [0.1 0.1 0.1 0.5]} ; 
end

sigmat = ll.sig.B.more.longest ; 
h = imsc_grid_comm(empmat.*sigmat ,remap_labs,extras{:}) ; 
h.AlphaData = sigmat ~= 0 ; 
xticks('')
axis square

if sdx ~= 1
    yticks('') ; 
end

xlabel('sig. more than surrogate')
title('ultra long')

nexttile()

sigmat = ll.sig.B.less.longest ; 
h = imsc_grid_comm(empmat.*sigmat ,remap_labs,extras{:}) ; 
h.AlphaData = sigmat ~= 0 ; 
xticks('')
axis square

xlabel('sig. less than surrogate')

if sdx ~= 1
    yticks('') ; 
end

%%

out_figdir = [ './reports/figures/supp/' ]
mkdir(out_figdir)
filename = [out_figdir '/surr_test_edgecounts.pdf' ] ; 
print(filename,'-dpdf','-vector')
close(gcf)

%%
% 
% filename = [DD.PROC '/surrogate3_' OUTSTR '_' , num2str(SPK_THR) , '_spk.mat'] ; 
% ll = load(filename) ; 