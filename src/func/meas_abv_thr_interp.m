function abvthr = meas_abv_thr_interp(ints,thr)

[~,~,~,aa] = get_contact_times(ints>thr) ; 

if isnan(aa)
    abvthr = nan ; 
    return 
end

ntp = length(ints) ;
nabv = size(aa,1) ; 
abvthr=nan(nabv,1) ; 

for idx = 1:nabv

    riseind = aa(idx,1) ; 
    fallind = aa(idx,2) ; 

    if riseind==1 
        riseadd=0 ; 
    else
        xx = [ riseind-1 riseind] ; 
        yy = [ ints(riseind-1) ints(riseind) ] ;
        riseadd = riseind-get_cross(xx,yy,thr); 
    end

    if fallind==ntp
        falladd=0 ; 
    else
        xx = [ fallind fallind+1] ; 
        yy = [ ints(fallind) ints(fallind+1) ] ;
        falladd = get_cross(xx,yy,thr) - fallind ; 
    end
     
    abvthr(idx) = fallind-riseind + riseadd + falladd ; 

end

end

function x0 = get_cross(xcross,ycross,thr) 
    % slope (rise/run)
    m = (ycross(2) - ycross(1)) / (xcross(2) - xcross(1) ) ;
    % y intersect
    b = ycross(1) - (m*xcross(1)) ;
    % find point where thr is crossed
    x0 = (thr-b)/m; 
end