function [bmat_ets] = get_blocky_ets(ets,ca)

% if size(indat,1) ~= length(ca) | size(indat,2) ~= length(ca)
%     error('wrong sizes')
% end

uniqca = unique(ca) ; 

%%

[ntp,nedges] = size(ets) ; 
nnodes = sqrt(2* nedges+0.25)-0.5+1 ; 

[u,v] = find(triu(ones(nnodes),1)) ; 

ca_u = ca(u) ;
ca_v = ca(v) ; 

%% run it

bmat_ets = zeros(length(uniqca),length(uniqca),ntp) ; 

for idx = 1:length(uniqca)

    inds_d = ca_u==idx & ca_v == idx ; 
    bmat_ets(idx,idx,:) = mean(ets(:,inds_d),2) ; 
    
    for jdx = idx+1:length(uniqca)
        
        inds_d2 = ca_u==jdx & ca_v == jdx ;

        inds_od = (ca_u == idx | ca_u == jdx) & (ca_v == idx | ca_v == jdx) ; 

        bmat_ets(idx,jdx,:) = mean(ets(:, ...
            logical(inds_od-inds_d-inds_d2) ),2) ; 
        bmat_ets(jdx,idx,:) = bmat_ets(idx,jdx,:) ; 
    
    end
end
%%

% bmat = cell2mat(...
%         arrayfun(@(i_) ...
%             arrayfun(@(j_) mean(ets(:,(ca_u==i_ & ca_v==j_)),2) , uniqca(:)', 'UniformOutput',false), ...
%         uniqca(:)',  'UniformOutput',false)' ...
%     ) ; 
