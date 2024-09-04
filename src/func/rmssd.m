function outdiff = rmssd(indat,k,difftype,cumu)
% root mean sucessive squared differences for each column
% 
% k defines order or lag of diff, based on difftype
% for difftype lag, have option to compute cumulative lag diff
%

if nargin < 2
    k = 1 ; 
end

if nargin < 3
    difftype = 'discrete' ; 
end

if nargin < 4
    cumu = false ; 
end

n = size(indat,1) ; 

if k == 1
    outdiff = (sum(diff(indat).^2)/(n-1)).^0.5 ; 

elseif strcmp(difftype,'discrete')
    outdiff = (sum(diff(indat,k).^2)/(n-k)).^0.5 ; 

elseif strcmp(difftype,'lag')
    if ~cumu % at just k
        outdiff = arrayfun(@(i_) (sum(diag(indat(:,i_) - indat(:,i_)',k).^2)/(n-1)).^0.5, ...
            1:size(indat,2) ) ; 
    else % mean for k and below
        outdiff = mean(cell2mat(arrayfun(@(j_) arrayfun(@(i_) (sum(diag(indat(:,i_) - indat(:,i_)',j_).^2)/(n-1)).^0.5, ...
            1:size(indat,2) ), 1:k,'UniformOutput',false)')) ; 
    end
else
    error('did not compute diff')
end

