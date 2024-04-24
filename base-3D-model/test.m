% Define the given vector and a point through which it passes
u = [2, 3]; % Example vector
A = [1, 1]; % Example point

% Find any other point on the line defined by the vector
B = A + u;

% Calculate the displacement vector
d = B - A;

% Find a vector perpendicular to d
p = [-d(2), d(1)];

% Translate the vector back to the desired position
B_new = B + p;

other = A - B_new;
disp(other)
% Plot the original vector and the perpendicular vector
figure;
plot([A(1), B(1)], [A(2), B(2)], 'b', 'LineWidth', 2); % Original vector
hold on;
plot([B(1), B_new(1)], [B(2), B_new(2)], 'r', 'LineWidth', 2); % Perpendicular vector
plot(A(1), A(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Point A
plot(B(1), B(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Point B
plot(B_new(1), B_new(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % New point B
plot([B_new(1), other(1)], [B_new(2), other(2)], 'g', 'LineWidth', 2)
grid on;
axis equal;
legend('Original Vector', 'Perpendicular Vector', 'Point A', 'Point B', 'New Point B');