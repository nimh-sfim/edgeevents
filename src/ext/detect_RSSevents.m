function  [pk_inds,pk_amp,numpk,pval,pcnt] = detect_RSSevents(time_series,Nreps,offsets,pthr)

Nchannels = size(time_series,2);
Ntp = size(time_series,1);

if nargin < 2
    Nreps = 1000 ;
end

if nargin < 3
    offsets = 1:Ntp ;
end

if nargin < 4
    pthr = 0.05 ; 
end

% compute ets and rss
ets = fcn_ets(time_series);
rssts = sum(ets.^2,2).^0.5;

% circshift null
pcnt = zeros(Ntp,1);
for r=1:Nreps % make the null distribution with this link
    disp(r)
    tsr = zeros(Ntp,Nchannels);
    for n=1:Nchannels
        % circshift function moves the timeseries left/right and shifts it
        % around to the front
        tsr(:,n) = circshift(time_series(:,n),offsets(randi(length(offsets))));
    end
    etsr = fcn_ets(tsr); % get null edge time series
    rsstsr = sum(etsr.^2,2).^0.5; % get rss of null edge time series
    pcnt = pcnt+(rssts>max(rsstsr));        % '1' if rssts>all rsstsr values
end

% pval
pval = 1-pcnt./Nreps;

% determine peaks as intersection of 'findpeaks' and 'pvals'
% findpeaks
[~, fp_ts] = findpeaks(rssts);
% intersection
pk_inds = intersect(find(pval<pthr),fp_ts);  
pk_amp = rssts(pk_inds);
numpk = length(pk_inds);

end

function [a,aa] = fcn_ets(ts)

[~,n] = size(ts);               % number samples/nodes
z = zscore(ts);                 % z-score
[u,v] = find(triu(ones(n),1));  % get edges
a = z(:,u).*z(:,v);             % edge ts products
if nargout>1
    aa = sqrt(z(:,u).^2 + z(:,v).^2); 
end
end