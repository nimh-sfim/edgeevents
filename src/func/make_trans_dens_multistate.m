function [ trans_mat ] = make_trans_dens_multistate(state_mat,tspan) 

if nargin < 2
    tspan = 2 ;
end

if tspan < 2
    error('tspan must be at least 2')
end

[n_states,ntp] = size(state_mat) ; 

trans_mat = zeros(n_states,n_states,n_states) ; 

for idx = 1:n_states

    curr_st = idx ; 
    curr_st_trans = zeros(n_states,n_states) ; 
    for tdx = 1:ntp-tspan+1
        if state_mat(curr_st,tdx) % if the curr state is active, look at the transitions
            cc = state_mat(:,(tdx+1):(tdx+tspan-1)) ; 
            curr_st_trans = ( (cc*cc')./ (tspan-1) ) + curr_st_trans ;  
        end
    end
    
    curr_st_trans = curr_st_trans ./ sum(state_mat(curr_st,:)) ; 
    trans_mat(:,:,idx) = curr_st_trans ; 

end
