function [ dens_mat ] = make_trans_dens(state_mat,tspan,countcurrent,summaryfunc) 

if nargin < 2
    tspan = 2 ; % tspan includes curr tp + how many tp in future
                % so tspan 2 == curr tp + 1 tp in future
end

if tspan < 1
    error('tspan must be at least 1')
end

if nargin < 3
    countcurrent = 0 ; 
end

if nargin < 4
    summaryfunc = @mean ; 
end

if ~countcurrent && (tspan < 2) 
    error('if not counting current, tspan needs to be at least 2')
end

%% run it

[n_states,ntp] = size(state_mat) ; 

dens_mat = zeros(n_states) ; 

for idx = 1:n_states

    curr_st = idx ; 
    curr_st_trans = zeros(n_states,1) ; 
    for tdx = 1:ntp-tspan+1 % loop through time points
        if state_mat(curr_st,tdx) > 0 % if the curr state is active, look at the states

            if countcurrent
                inds = (tdx):(tdx+tspan-1) ;
            else
                inds = (tdx+1):(tdx+tspan-1) ;
            end

            curr_st_trans = summaryfunc(state_mat(:,inds),2) + curr_st_trans; 
        end
    end
    
    % % normalize relative to how many each state spikes
    % curr_st_trans = curr_st_trans ./ sum(state_mat(curr_st,:)) ;                                                                 
    dens_mat(idx,:) = curr_st_trans ; 

end


