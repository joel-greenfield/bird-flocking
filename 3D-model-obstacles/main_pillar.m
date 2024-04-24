%% Agent-based model of flocking birds
%% Cylindrical pillar boundary

tic

% Parameters
n = 200; % Number of birds in flock
p_f2l_0 = 2e-4; % Initial prob. of follower becoming leader (geometric distribution)
d_per = 20; % Persistence distance (max. distance to nearest neighbor as leader)
t_per = 700; % Persistence time (max. time as leader)
t_ref = 800; % Refractory time (min. time as follower before becoming leader)
t_del = 0.1; % Delay time (reaction time in velocity delay diff. eq.)
C_ali = 3; % Alignment coefficient
C_att = 0.01; % Attraction coefficient
C_rep = 2.5; % Repulsion coefficinet
ep = 1e-6; % Small number to avoid degenerate repulsion forces
M = 7; % Number of nearest neighbors affecting birds


% Initial conditions
x_0 = 200*rand(n,3); % Position (unif. distributed in cube w/ side length 200)
v_0 = zeros(n,3); % Velocity (all zero)

% Iterate model
t_0 = 0; % Initial time
dt = 0.1; % Time increment
T = 2e3; % Final time

t = t_0; % Time

% Initialized vectors
status = false(n,1); % Status (initialized as all followers; 0 = follower, 1 = leader)
status_h = zeros(n,T/dt); status(:,1) = status; % History of each bird's status
p_f2l = p_f2l_0*ones(n,1); % Prob. of follower becoming leader (geometric distribution; initialized as all initial prob.)
t_status = t_ref*ones(n,1); % Time as current status (initialized as all followers after refractory time)
x = x_0; % Position
v = v_0; % Velocity
x_h = zeros(n,3,T/dt); x_h(:,:,1) = x_0; % History of each bird's position
v_h = zeros(n,3,T/dt); v_h(:,:,1) = v_0; % History of each bird's velocity

% % %
n_f = zeros(T/dt,1); % Number of followers
% % %

while t < T
    f = (status == 0); % Follower status (0 = leader, 1 = follower)
    l = (status == 1); % Leader status (0 = follower, 1 = leader)
    
    % % %
    n_f(round(t/dt)+1) = sum(f);
    % % %

    %% Followers become leaders
    f2l_1 = (rand(n,1) <= p_f2l); % Stochastic process (geometric distribution)
    f2l_2 = (t_status >= t_ref); % Time as follower exceeds refractory time
    f2l = f & f2l_1 & f2l_2; % Followers becoming leaders

    status(f2l) = 1; % Update status of new leaders
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

    status(l2f) = 0; % Update status of new followers
    p_f2l(l2f) = p_f2l_0; % Prob. of follower becoming leader set to initial prob.
    t_status(l2f) = 0; % Reset time as current status of new followers

    %% Keep birds within boundary
    pillar_centers = [[ 0      1000 0]; ...
                      [ 1000   0    0]; ...
                      [ 0     -1000 0]; ...
                      [-1000   0    0]];
    x = pillar_boundary(x,pillar_centers,200);

    %% Update variables
    p_f2l(t_status >= t_ref) = p_f2l(t_status >= t_ref)*(1-p_f2l_0); % Increment prob. of follower becoming leader for followers beyond refractory time (geometric distribution)
    t = t+dt; % Increment time
    t_status = t_status+dt; % Increment time as status
    status_h(:,round(t/dt)+1) = status;
    x_h(:,:,round(t/dt)+1) = x;
    v_h(:,:,round(t/dt)+1) = v;

    %% Update position & velocity
    v = v+velocity(f,l,x,x_h,v_h,t,t_del,dt,C_ali,C_att,C_rep,ep,M,status)*dt;
    x = x+v*dt;
end

toc

close all

animate_birds(x_h,dt,'testvideo')

figure
plot(n-n_f)

plot_centers(x_h)

function x = pillar_boundary(x,pillar_centers,pillar_radius)
    for i = 1:size(pillar_centers,1)
        pillar_center = pillar_centers(i,:);
        x1 = x-pillar_center;
        [x_t,x_r,x_z] = cart2pol(x1(:,1),x1(:,2),x1(:,3));
        x2 = [x_t x_r x_z];
        condition = x_r<pillar_radius;
    
        x2(condition,:) = [x_t(condition) pillar_radius*ones(size(x_r(condition))) x_z(condition)];
    
        [xnew1,xnew2,xnew3] = pol2cart(x2(:,1),x2(:,2),x2(:,3));
        x = [xnew1 xnew2 xnew3]+pillar_center;
    end
end