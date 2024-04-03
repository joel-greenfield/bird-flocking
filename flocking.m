%% Agent-based model of flocking birds

% Parameters
n = 200; % Number of birds in flock
p_f2l_0 = 2e-4; % Initial prob. of follower becoming leader (geometric distribution)
d_per = 20; % Persistence distance (max. distance to nearest neighbor as leader)
t_per = 700; % Persistence time (max. time as leader)
t_ref = 800; % Refractory time (min. time as follower before becoming leader)
t_del = 0.1; % Delay time (reaction time in velocity delay diff. eq.)

% Initial conditions
x_0 = 100*(2*rand(n,3)-1); % Position (unif. distributed in cube w/ side length 200)
v_0 = zeros(n,3); % Velocity (all zero)

% Initialized vectors
status_c = false(n,1); % Status (initialized as all followers; 0 = follower, 1 = leader)
status_h = zeros(n,T/dt); % History of each bird's status
p_f2l = p_f2l_0*ones(n,1); % Prob. of follower becoming leader (geometric distribution; initialized as all initial prob.)
t_status = t_ref*ones(n,1); % Time as current status (initialized as all followers after refractory time)
x = x_0; % Position
v = v_0; % Velocity

% Iterate model
t_0 = 0; % Initial time
dt = 0.1; % Time increment
T = 3e3; % Final time

t = t_0; % Time

% % %
n_f = zeros(T/dt,1); % Number of followers
% % %

while t <= T
    f = (status_c == 0); % Follower status (0 = leader, 1 = follower)
    l = (status_c == 1); % Leader status (0 = follower, 1 = leader)
    
    % % %
    n_f(round(t/dt)+1) = sum(f);
    % % %

    %% Followers become leaders
    f2l_1 = (rand(n,1) <= p_f2l); % Stochastic process (geometric distribution)
    f2l_2 = (t_status >= t_ref); % Time as follower exceeds refractory time
    f2l = f & f2l_1 & f2l_2; % Followers becoming leaders

    status_c(f2l) = 1; % Update status of new leaders
    p_f2l(f2l) = 0; % Prob. of follower becoming leader set to zero for leaders
    t_status(f2l) = 0; % Reset time as current status of new leaders

    %% Leaders become followers
    l2f_1 = (t_status >= t_per); % Time as leader exceeds persistence time
    
    % Distance to nearest neighbor exceeds persistence distance
    l2f_2 = false(n,1);
    for i = 1:n
        d_neighbors = sqrt(sum((x-x(i,:)).^2,2));
        d_neighbors(i) = Inf;
        if min(d_neighbors) >= d_per
            l2f_2(i) = 1;
        end
    end

    l2f = l & (l2f_1 | l2f_2); % Leaders becoming followers

    status_c(l2f) = 0; % Update status of new followers
    p_f2l(l2f) = p_f2l_0; % Prob. of follower becoming leader set to initial prob.
    t_status(l2f) = 0; % Reset time as current status of new followers

    %% Update position & velocity
    v = v+dvdt(status_h)*dt;
    x = x+v*dt;

    %% Update variables
    p_f2l(t_status >= t_ref) = p_f2l(t_status >= t_ref)*(1-p_f2l_0); % Increment prob. of follower becoming leader for followers beyond refractory time (geometric distribution)
    status_h(:,round(t/dt)+1) = status_c;
    t = t+dt; % Increment time
    t_status = t_status+dt; % Increment time as status
end

% % %
close all
plot(n-n_f)
% % %

function dvdt = dvdt(status_h) %INCOMPLETE
    dvdt = 0;
end