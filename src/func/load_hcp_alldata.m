function [outStr] = load_hcp_alldata(dirStr,parc,sublist,globstr)
% dirStr = './data/hcp7t_test_LP/100610/' ;
% parc = yan200 ;

if nargin < 4
    globstr = '*MRI*' ;
end

% dd = dir(strcat( dirStr , '/*/' , globstr , '/*' , parc , '*nii.gz' ) ) ; 
% if isempty(dd)
%     error('did not match any data')
% end

outStr = struct();

if ischar(sublist) || isstring(sublist)
    aa = cellstr(num2str(load(sublist))) ; 
else
    aa = sublist ;
end

for idx = 1:length(aa) 

    dd = dir(strcat( dirStr , '/', aa{idx} , '/', globstr , '/*' , parc , '*nii.gz' ) ) ; 
    if length(dd) > 1
        error('globbed more than 1 match per person')
    end
    [~,b] = fileparts(dd.folder) ;
%     outStr.(matlab.lang.makeValidName(b)).ts = ...
%         squeeze(niftiread([ dd(idx).folder '/' dd(idx).name ])) ;
    outStr(idx).ts = squeeze(niftiread([ dd.folder '/' dd.name ])) ;
    outStr(idx).name = b ;
    outStr(idx).sub = aa{idx} ; 
end
