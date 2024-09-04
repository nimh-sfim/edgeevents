function [state_mat] = get_ets_states(ets,ca,ondiagthr)

if nargin < 3
    ondiagthr = 0.2 ; 
end

%% check the ca

if length(unique(ca)) ~= max(ca)
    error('please condition ca to have numbes 1:Nstates, no gaps')
end

%% setup

[ntp,nedges] = size(ets) ; 
nnodes = sqrt(2* nedges+0.25)-0.5+1 ; 

if length(ca) ~= nnodes
    error('ca wrong size')
end

%% run it

od_nan = @(x_) x_+single(diag(nan(size(x_,1),1))) ; 

% no loops
% first block mat so don't have to double make blocks
ca_blocks =  get_blocky_ets(ets,ca) ; 

% replacing 0's makes the following operations easier
ca_blocks(ca_blocks==0) = nan ; 

% get states trace, with on-diag threshold
if strcmp(ondiagthr,'assort') 
    % make sure the on diagonal is greater/equal to than the mean off diag
    state_mat = cell2mat(arrayfun(@(i_) ...
        sum(ca_blocks(:,:,i_),2,'omitnan').*single(diag(ca_blocks(:,:,i_))>=...
            mean(od_nan(ca_blocks(:,:,i_)),2,'omitnan')),1:ntp,'UniformOutput',false)) ;
elseif strcmp(ondiagthr,'assortweak') 
    % make sure the on diagonal is greater/equal to than the minimum off diag
    state_mat = cell2mat(arrayfun(@(i_) ...
        sum(ca_blocks(:,:,i_),2,'omitnan').*single(diag(ca_blocks(:,:,i_))>=...
            min(od_nan(ca_blocks(:,:,i_)),[],2)),1:ntp,'UniformOutput',false)) ;
elseif strcmp(ondiagthr,'assortstrong') 
    % make sure the on diagonal is greater/equal to than the max off diag
    state_mat = cell2mat(arrayfun(@(i_) ...
        sum(ca_blocks(:,:,i_),2,'omitnan').*single(diag(ca_blocks(:,:,i_))>=...
            max(od_nan(ca_blocks(:,:,i_)),[],2)),1:ntp,'UniformOutput',false)) ;
elseif ondiagthr >= 0 
    state_mat = cell2mat(arrayfun(@(i_) ...
        sum(ca_blocks(:,:,i_),2,'omitnan').*single(diag(ca_blocks(:,:,i_))>ondiagthr),1:ntp,'UniformOutput',false)) ;
else
    error('dont know what to do with ondiagthr')
end