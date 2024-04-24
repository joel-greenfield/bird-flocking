%% Calculate dVdt for one time step
function dvdt = velocity(f,l,x,x_h,v_h,t,t_del,dt,C_ali,C_att,C_rep,ep,M)
    temp = zeros(200, 3);
    if t>=t_del
        x_del = x_h(:,:,round((t-t_del)/dt)+1);
        v_del = v_h(:,:,round((t-t_del)/dt)+1);
    
        %% Find each bird's M nearest neighbors
        dist = zeros(length(f)); % dist(i,j) = distance between birds i & j
        for i = 1:length(f)
            dist(i,:) = sqrt(sum((x_del-x_del(i,:)).^2,2));
        end
        dist(dist==0) = Inf;
        [~,dist_index] = sort(dist,2);
        n_n = dist_index(:,1:M); % Each bird's M nearest neighbors
        
        % calculate dVdt for followers
        f_v = alignment(f,l,x_del,v_del,C_ali,C_att,C_rep,ep,n_n,M) ...
            +attraction(f,l,x_del,v_del,C_ali,C_att,C_rep,ep,n_n,M) ...
            +repulsion(f,l,x_del,v_del,C_ali,C_att,C_rep,ep,n_n,M);

        % calculate dVdt for leaders
        l_v = repulsion(f,l,x_del,v_del,C_ali,C_att,C_rep,ep,n_n,M);
        
        % set dvdt for all agents
        temp(f, :) = f_v(f,:); 
        temp(l, :) = l_v(l,:);
        dvdt = temp;


        
    else % intialize velocity to zero if time is less than time delay
        dvdt = zeros(length(f),3);
    end
end