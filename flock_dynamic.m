clear all;
global time_step;
global birds; 

bird = struct('X', double.empty(0, 3), 'V', double.empty(0, 3), 's_k', 0, 'leader_time', 0);
% X and V are Nx2 arrays where X_ij is the position and velocity (3 dim)
% at time t = i - 1

N = 20;    
birds = repmat(bird(), N, 1);

%birds = createArray(200, 1, bird);
init_Xx = randperm(100, 20); % 200 birds within 20 ft? apart
init_Xy = randperm(100, 20); 
init_Xz = randperm(100, 20); 

for i=1:size(birds, 1)
    % initialize bird position at t=0 to be random
    birds(i).X(1,:) = [init_Xx(i), init_Xy(i), init_Xz(i)];
    birds(i).V(1,:) = [0, 0, 0]; % initial velocity 
end
time_step = 0; % start 0 time step
for t=2:35
    time_step = time_step + 1;
    flock_dynamics(7, t)
end

fig = gcf;
for i=1:time_step
    disp(i)
    clf(fig)
    ax = fig.CurrentAxes;
    for j=1:size(birds,1)
        if birds(j).s_k == 1
            scatter3(birds(j).X(i, 1), birds(j).X(i, 2), ...
                birds(j).X(i, 3), 'filled', 'MarkerFaceColor', 'red')
        else
            scatter3(birds(j).X(i, 1), birds(j).X(i, 2), ...
                birds(j).X(i, 3), 'filled', 'MarkerFaceColor', 'black')
        end
        hold on
        xlim([0 2000])
        ylim([0 2000])
        zlim([0 2000])
    end
    hold off
    pause(0.1)
end

% try to solve odes with constant nearest neighbors for 0 to t. This time
% will be seperate from the time used to update the birds position: ie
% call ode45 on this for each time step and we will get a vector for the
% birds X and V. Then what values do we update X and V with?
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
    delta = 0;

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
    delta = 0;

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



function flock_dynamics(M, time)
    global time_step;
    global birds;
    % C_rep = 2.5;
    % C_ali = 3;
    % C_att = 0.01;
    delta = 0;
    p = 3; % persistance time (original 700)
    d = 20; % persistance distance
    tau = 800;
    p_transition = .0002;
    
    M_count = 0;
    
    for i=1:size(birds, 1)
        random_num = randi(1000, 'double');
        if random_num <= 5
            birds(i).s_k = 1;
            disp("BIRD " + i + " IS LEADER")
        end
        if birds(i).s_k == 1
            birds(i).leader_time = birds(i).leader_time + 1;
            if birds(i).leader_time > p % bird has been leader for longer time than persistance time
                birds(i).s_k = 0;    % reset to follower
                birds(i).leader_time = 0;
            end
        end
        dist = zeros(size(birds, 1), 1);
        NN_idx = zeros(M, 1); % indices of the M nearest neighbors to bird i
        for j=1:size(birds, 1) 
            if i ~= j
                if time_step == 1 % first time step, compare everyone's inital positions
                    dist(j) = norm(birds(i).X(time_step,:) - birds(j).X(time_step,:));
                else    % some birds have been updated, want to compare previous time steps
                    dist(j) = norm(birds(i).X(time_step - 1,:) - birds(j).X(time_step - 1,:));
                end
            else
                dist(j) = Inf; % distance to self
            end
            [~, idx] = sort(dist);
            NN_idx = idx(1:M, :); % splice to get the M nearest neighbors' indices
            if dist(NN_idx(1)) > d % distance to nearest neighbor is greater than persistance dist
                birds(i).s_k = 0;
                birds(i).leader_time = 0;
            end
        end
        delta = 0;
        % make array of this birds positions and velocities at time t - delta

        y0 = [birds(i).X(time_step - delta,:) birds(i).V(time_step - delta,:)];
        tspan = [0 .1];
        [t, y] = ode45(@(t,y) bird_ODEs2(t,y,NN_idx, birds(i)),tspan,y0);
        % update X and V for this time step
        birds(i).X(time_step + 1, :) = y(end, 1:3);
        birds(i).V(time_step + 1, :) = y(end, 4:6);
        %disp(y)
        %disp(birds(i).X)    %pass by reference?? value of birds(i) changes after exiting this loop
    end
end