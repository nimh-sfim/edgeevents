function [ dist_mat ] = make_tinda(state_mat) 

%% run it

[n_states,ntp] = size(state_mat) ; 

dist_mat = zeros(n_states) ; 

for idx = 1:n_states

    firsthalf_dens = zeros(n_states,1) ; 
    secndhalf_dens = zeros(n_states,1) ;

    curr_st = idx ; 
    count_intervals = 0 ; 
    for tdx = 1:ntp % loop through time points
        if state_mat(curr_st,tdx) > 0 % if the curr state is active, look at the states

            % find next spike
            nxtspk =  tdx + find(state_mat(curr_st,tdx+1:end),1,"first") ; 
            if isempty(nxtspk)
                % need to have another spike to compute
                continue
            end

            tspan = [tdx+1 nxtspk-1 ] ;
            if range(tspan) < 4 
                % need at least 4 to differentiate
                % disp('not 4')
                continue
            end

            midind = ( range(tspan)/2 ) + tspan(1)  ;  

            firsthalf_ind = tspan(1):floor(midind) ; 
            secndhalf_ind = ceil(midind):tspan(2) ;

            firsthalf_dat = state_mat(:,firsthalf_ind) ;
            secndhalf_dat = state_mat(:,secndhalf_ind) ;

            % % if even number, divide appropriate column by 2
            % if mod(midind,1)==1
            %     firsthalf_dat(:,end) = firsthalf_dat(:,end) ./ 2 ; 
            %     secndhalf_dat(:,1) = secndhalf_dat(:,1) ./ 2 ; 
            % end

            firsthalf_dens = firsthalf_dens + mean(firsthalf_dat,2) ; 
            secndhalf_dens = secndhalf_dens + mean(secndhalf_dat,2) ; 

            count_intervals = count_intervals + 1  ; 

            % dist_mat(curr_st,:) = dist_mat(curr_st,:) + firsthalf_dens(:)' ; % make horizontal
            % dist_mat(:,curr_st) = dist_mat(:,curr_st) + secndhalf_dens ; 

            % dist_mat(curr_st,:) = dist_mat(curr_st,:) + ...
            %     single((firsthalf_dens(:)' - secndhalf_dens(:)')>0) ; 

        end
    end

    dist_mat(curr_st,:) = ( firsthalf_dens(:)' - secndhalf_dens(:)' ) ./ count_intervals ; 
    
end
