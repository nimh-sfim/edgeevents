function [bmat] = get_blocky(indat,ca,summaryfunc)

if size(indat,1) ~= length(ca) | size(indat,2) ~= length(ca)
    error('wrong sizes')
end

if nargin < 3
    summaryfunc = @mean ; 
end

uniqca = unique(ca) ; 

bmat = cell2mat(...
        arrayfun(@(i_) ...
            arrayfun(@(j_) summaryfunc(unroll_helper(indat(ca==i_,ca==j_),i_,j_),'all') , uniqca(:)' ), ...
        uniqca(:)',  'UniformOutput',false)' ...
    ) ; 

% N = length(uniqca) ; 
% bmat = zeros(N) ; 
% for idx = 1:N
%     for jdx = idx:N
%         i_ = uniqca(idx) ; 
%         j_ = uniqca(jdx) ; 
%         aaa(idx,jdx) = sumfunc(indat(ca==i_,ca==j_),'all') ; 
%     end
% end
% bmat = triu(aaa)+triu(aaa,1)' ; 

end

% make a lil helper function, to handle self-loops that we don't want to
% look at!! 
function omat = unroll_helper(inmat,i1,i2)

    if i1~=i2
        omat = inmat(:) ;
    else
        mask = logical(triu(ones(size(inmat)),1));
        omat = inmat(mask) ;
    end
end