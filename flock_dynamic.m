clear all;
global time_step;
global birds; 

%% bird struct

% X and V: Nx3 arrays where X_i is position at time step i and 
%          V_i is velocity at time step i
% s_k: vector of length time_step, where s_k[i] is the bird's status at
%      time step i. (follower = 0, leader = 1)
% t_leader: number of time steps bird has been a leader
% t_follower: number of time steps bird has been a follower
% p_leader: bird's probability of becoming a leader
bird = struct('X', double.empty(0, 3), 'V', double.empty(0, 3), 's_k', uint8.empty(0, 1), ...
    't_leader', 0, 't_follower', 0, 'p_leader', 2e-4);

%% model flock dynamics for N birds

N = 20;    
birds = repmat(bird(), N, 1); % create vector of N birds

init_pos = 100*(2*rand(N,3)-1) + 500; % Nx3 vector of random positions

init_Xx = randperm(100, 20); % 200 birds within 20 ft? apart
init_Xy = randperm(100, 20); 
init_Xz = randperm(100, 20); 

% initialize bird structs
for i=1:size(birds, 1)
    birds(i).X(1,:) = init_pos(i, :);
    birds(i).V(1,:) = [0, 0, 0]; % initial velocity 
    birds(i).s_k(1,:) = 0;
end

T = 1000; % 100 seconds
% loop for T time steps (each cooresponding to dt = 0.1)
time_step = 2; % start 2nd time step (time step 1 is initial conditions)
for t=2:T
    flock_dynamics(4)
    time_step = time_step + 1;
end

%% plot solutions for each time step
fig = gcf;
for i=1:(time_step - 1) / 5
    disp(i)
    clf(fig)
    ax = fig.CurrentAxes;
    for j=1:size(birds,1)
        if birds(j).s_k(i) == 1
            scatter3(birds(j).X(i*5, 1), birds(j).X(i*5, 2), ...
                birds(j).X(i*5, 3), 'filled', 'MarkerFaceColor', 'red')
        else
            scatter3(birds(j).X(i*5, 1), birds(j).X(i*5, 2), ...
                birds(j).X(i*5, 3), 'filled', 'MarkerFaceColor', 'black')
        end
        hold on
        xlim([0 1000])
        ylim([0 1000])
        zlim([0 1000])
    end
    hold off
    pause(0.0001)
end

%% ODE func

% try to solve odes with constant nearest neighbors for each time step. This time
% will be seperate from the time used to update the birds position: ie
% call ode45 on this for each time step and we will get a vector for the
% birds X and V at t = timestep * dt
function dydt = bird_ODEs2(t, y, NN_idx, curr_bird)
    global time_step;
    global birds;
    X = y(1:3)'; % should be this birds postition at time t - delta
    %curr_bird.X(t - delta,:)
    V = y(4:6)'; % should be this birds velocity
    %curr_bird.V(t - delta,:)

    % parameters
    C_rep = 2.5;
    C_ali = 3;
    C_att = 0.01;
    epsilon = 0.01;
    delta = 1;

    % make array of birds positions and velocities at time t - delta
    M = size(NN_idx, 1);
    neighbor_pos = zeros(M, 3);
    neighbor_vel = zeros(M, 3);
    for i=1:M
        neighbor_pos(i, :) = birds(NN_idx(i)).X(time_step - delta, :);
        neighbor_vel(i, :) = birds(NN_idx(i)).V(time_step - delta, :);
    end

    % numerator = neighbor_pos - X;
    % denom = norm(neighbor_pos - X)^2 + epsilon;
    % A_rep = -C_rep * sum(numerator / denom);
    denom = vecnorm(neighbor_pos, 2, 2).^2 + norm(X, 2)^2 - 2*(neighbor_pos*X');
    A_rep = -1*C_rep * sum((neighbor_pos - X) ./ ((vecnorm(neighbor_pos, 2, 2).^2 + norm(X, 2)^2 - 2.*(neighbor_pos*X')) + epsilon));
    if curr_bird.s_k(time_step) == 1 % if leader at time_step
        dVdt = A_rep;
    else
        A_ali = (C_ali/M) * sum(neighbor_vel - V);
        A_att = C_att * sum(neighbor_pos - X);
        dVdt = A_rep + A_ali + A_att;
    end
    dXdt = V;
    dydt = [dXdt dVdt]';
end


% i want to solve this ode one time step at a time so that I can update
% the current bird's position and velocity (outside of ode45)
function dydt = bird_ODEs(t, y, NN_idx, birds, curr_bird)
    global time_step;
    time_step = time_step + 1;
    X = y(1:3)'; % should be this birds postition at time t - delta
    %curr_bird.X(t - delta,:)
    V = y(4:6)'; % should be this birds velocity
    %curr_bird.V(t - delta,:)
    % parameters
    C_rep = 2.5;
    C_ali = 3;
    C_att = 0.01;
    epsilon = 0.01;
    delta = 1;

    % make array of birds positions and velocities at time t - delta
    M = size(NN_idx, 1);
    neighbor_pos = zeros(M, 3);
    neighbor_vel = zeros(M, 3);
    for i=1:M
        neighbor_pos(i, :) = birds(NN_idx(i)).X(time_step - delta, :);
        neighbor_vel(i, :) = birds(NN_idx(i)).V(time_step - delta, :);
    end

    numerator = neighbor_pos - X;
    denom = norm(neighbor_pos - X)^2 + epsilon;
    A_rep = -C_rep * sum(numerator / denom);
    if curr_bird.s_k == 1
        dVdt = A_rep;
    else
        A_ali = (C_ali/M) * sum(neighbor_vel - V);
        A_att = C_att * sum(neighbor_pos - V);
        dVdt = A_rep + A_ali + A_att;
    end
    dXdt = V;
    dydt = [dXdt dVdt]';
end



function flock_dynamics(M)
    global time_step;
    global birds;
    delta = 1; % # of time steps for delay (1 time step: dt = 0.1)
    p = 700; % persistance time (original 700)
    d = 200; % persistance distance
    refractory_time = 800;
    p_transition = 2e-4;
    
    
    for i=1:size(birds, 1)
        random_num = rand(); % generate random double btween 0 and 1

        % switch to leader at probability p_leader
        if random_num <= birds(i).p_leader 
            assert(birds(i).s_k(time_step - 1) == 0)
            birds(i).s_k(time_step) = 1;
            birds(i).p_leader = 0;
            birds(i).t_leader = birds(i).t_leader + 1;
            birds(i).t_follower = 0;
            disp("BIRD " + i + " IS LEADER")

        % bird was a leader on prev. time step
        elseif birds(i).s_k(time_step - 1) == 1 
            if birds(i).t_leader > p % bird has been leader longer than persistance time
                % reset to follower
                birds(i).s_k(time_step) = 0;
                birds(i).t_leader = 0;
                birds(i).t_follower = birds(i).t_follower + 1;
                birds(i).p_leader = p_transition; % reset probability to initial probability
            else
                birds(i).s_k(time_step) = 1; % stay a leader
                birds(i).t_leader = birds(i).t_leader + 1;
            end
            
        % was a follower on prev. time step and did not switch status this time step
        else 
             birds(i).s_k(time_step) = 0; % stay a follower
             birds(i).t_follower = birds(i).t_follower + 1;
             if birds(i).t_follower > refractory_time % bird has been follower for longer than refractory time
                birds(i).p_leader = birds(i).p_leader * p_transition; % decrease prob. of becoming leader
             end
        end


        distances = zeros(size(birds, 1), 1);
        NN_idx = zeros(M, 1); % indices of the M nearest neighbors to bird i
        for j=1:size(birds, 1) 
            if i ~= j
                distances(j) = norm(birds(i).X(time_step - 1,:) - birds(j).X(time_step - 1,:));
            else
                distances(j) = Inf; % distance to self
            end
            [~, idx] = sort(distances);
            NN_idx = idx(1:M, :); % splice to get the M nearest neighbors' indices
        end
        
        % distance to nearest neighbor is greater than persistance dist
        if distances(NN_idx(1)) > d & birds(i).s_k(time_step - 1) == 1
            % reset to follower   
            birds(i).s_k(time_step) = 0;
            birds(i).t_leader = 0;
            birds(i).t_follower = birds(i).t_follower + 1;
            birds(i).p_leader = p_transition;
        end
        % make array of this birds position and velocity at time t - delta
        y0 = [birds(i).X(time_step - delta,:) birds(i).V(time_step - delta,:)];
        tspan = [0 0.1]; % [0 dt]
        [t, y] = ode45(@(t,y) bird_ODEs2(t,y,NN_idx, birds(i)),tspan,y0);
        %y = dde23(@(t,y) bird_ODEs2(t, y, NN_idx, birds(i)), 0.1, @bird_hist, tspan);
        % update X and V for this time step
        birds(i).X(time_step, :) = y(end, 1:3);
        birds(i).V(time_step, :) = y(end, 4:6);
        %disp(y)
        %disp(birds(i).X)    %pass by reference?? value of birds(i) changes after exiting this loop
    end
end

function history = bird_hist(t)
    history = ones(20, 6);
end