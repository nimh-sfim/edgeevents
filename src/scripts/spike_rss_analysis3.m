%% start

clearvars
clc
close all

%% preamble load data

run('./config/config_hcp_sch200_1.m') 

SPK_THR = 2.25 ;

high_bin = 4 ; 
maxspk = 1200  ; 
lowmedhigh_edges = [ 1 2 high_bin maxspk ] ; 

spklen_names = {'short' 'inter' 'long'} ; 

% mask to only pickout cortical nodes
[u,v] = find(triu(ones(255),1));  
cortmask =  u<=finfo.nnodes & v<=finfo.nnodes  ;

subsets = {'subset1' 'subset2'} ; 
nedges = ( finfo.nnodes * (finfo.nnodes-1) ) / 2 ; 

edgetrim = 42 ;
timevec = 1:finfo.ntp ;  
wanttime = timevec((edgetrim/2):end-(edgetrim/2)-1) ;


%% get the empirical slopes and the kendall

for sdx = subsets

    for idx = 1:length(sublist.(sdx{1}))
    
        disp([ num2str(idx) '-' num2str(length(sublist.(sdx{1}))) ' ' sdx{1}])  ; 
    
        sind = find(cellfun(@(x_)strcmp(x_,sublist.(sdx{1})(idx)),sublist.all)) ; 
    
        filename = [DD.PROC '/' imglob '/' datStr(sind).sub '_' OUTSTR '_' , num2str(SPK_THR) , '_spike_len.mat'] ; 
        readdat = load(filename,'spike_len_mat') ; 
    
        ets = get_ets(datStr(sind).ts(wanttime,1:finfo.nnodes)) ; 
            
        dd = discretize( readdat.spike_len_mat(:,cortmask), lowmedhigh_edges) ; 
        dd(isnan(dd)) = 0 ;

        mm = cell(3,1) ; 
        for jdx = 1:3
            [~,oo] = count_spks(dd==jdx) ;  
            ss = sum(oo(wanttime,:),2) ;

            [~,ii] = spk_pks_sort(ss,'MinPeakProminence',std(ss)./2) ; 
            
            ll = min(50,length(ii)) ; 

            mm{jdx} = mksq(mean(ets(ii(1:ll),:))) ; 
               
        end

        [~,o1] = count_spks(dd==1) ; 
        [~,o2] = count_spks(dd==2) ; 
        [~,o3] = count_spks(dd==3) ; 

        

        sum(o1(wanttime,:),2) >= prctile(o1(wanttime,:),95)
        


    end
end