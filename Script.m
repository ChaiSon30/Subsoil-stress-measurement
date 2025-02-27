clear all
load("Pressurefunction.mat") % load pressure distribution prediction function

% Given parameters
z = 0.0435;  % Sinkage (m)
l = 0.3939;  % Contact patch length (m)
R = 0.785/2;  % Radius of tire (m)

% Define the system of equations
fun = @(theta) [ (cos(theta(2)) - cos(theta(1))) - (z/R);
                 (sin(theta(1)) + sin(theta(2))) - (l/R) ];

% Initial guess for theta_e and theta_r
theta0 = [pi/4, pi/4]; 

% Solve for theta_e and theta_r
options = optimset('Display', 'iter'); % Show iterations
theta_solution = fsolve(fun, theta0, options);

% Extract solutions
theta_e = theta_solution(1);
theta_r = theta_solution(2);

theta_a = rad2deg(theta_e);
theta_r = rad2deg(theta_r);

%% Assign depth values to contact patch

% Parameters
b = 0.104; % Semi-minor axis of superellipse (m)
n = 3; % Power for superellipse
rotation_angle = 0; % Rotation angle (degrees)
square_size = 0.01; % Length of smaller squares (m)
depth = 0.15; % Depth below the superellipse (m)
R = 0.785/2;
a = (R*sin(deg2rad(theta_a)) + R*sin(deg2rad(theta_r)))/2;
d = R-R*cos(deg2rad(theta_a));

% Generate grid for the superellipse's bounding box
x_range = -3*a:square_size:3*a;
y_range = -3*a:square_size:3*a;
[X, Y] = meshgrid(x_range, y_range);

% Rotate the grid points
theta = deg2rad(rotation_angle);
X_rot = X * cos(theta) - Y * sin(theta);
Y_rot = X * sin(theta) + Y * cos(theta);

% Superellipse equation
inside_superellipse = (abs(X_rot/a).^n + abs(Y_rot/b).^n) <= 1;
% Assign pressure to squares
pressure = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        if inside_superellipse(i, j)
            vertices_x(i,j) = X_rot(i,j);
            vertices_y(i,j) = Y_rot(i,j);

            if (-a <= sqrt(X_rot(i,j)^2 + Y_rot(i,j)^2)) && (sqrt(X_rot(i,j)^2 + Y_rot(i,j)^2)) <= (R*sin(deg2rad(theta_a))-a)
                vertices_z(i,j) = (R - R*cos(asin((R*sin(deg2rad(theta_a))-a-X_rot(i,j))/R)))-d;
            elseif (sqrt(X_rot(i,j)^2 + Y_rot(i,j)^2) > (R*sin(deg2rad(theta_a))-a)) && (sqrt(X_rot(i,j)^2 + Y_rot(i,j)^2) <= a)
                vertices_z(i,j) = (R-R*cos(asin((X_rot(i,j)-R*sin(deg2rad(theta_a))+a)/R)))-d;
            elseif sqrt(X_rot(i,j)^2 + Y_rot(i,j)^2) > a
                vertices_z(i,j) = -(R - R*cos(deg2rad(theta_a)))+(R - R*cos(deg2rad(theta_r)));
            end
            yd(i,j) = (R - R * cos (atan(vertices_y(i,j)/R)));
            vertices_z(i,j) = vertices_z(i,j) + yd(i,j);
        end
    end
end

zero_indices = (vertices_x == 0 & vertices_y == 0 & vertices_z == 0);
vertices_z(zero_indices) = NaN;

figure;
surf(vertices_x, vertices_y, vertices_z, 'FaceColor', 'interp', 'EdgeColor', 'none');

% Use colormap to indicate depth
colormap("parula"); 
cb = colorbar; 
ylabel(cb, "Depth (m)", 'Rotation',270)
caxis([min(vertices_z(:)), max(vertices_z(:))]); % Scale color range to z-values

%% Plot contact patch

% Formatting the plot
shading interp; % Smooth color shading
grid on;
view(3); % Set 3D view
xlabel('X co-ordinate (m)');
ylabel('Y co-ordinate (m)');
zlabel('Z co-ordinate (m)');
axis equal
hold on
theta1 = linspace(0, 2*pi, 100);
ellipse_x = ((abs(cos(theta1)/a)).^3 + (abs(sin(theta1)/b)).^3).^(-1/3) .* cos(theta1);
ellipse_y = ((abs(cos(theta1)/a)).^3 + (abs(sin(theta1)/b)).^3).^(-1/3) .* sin(theta1);
ellipse_z = zeros(length(theta1)); % All points on z=0 plane
plot3(ellipse_x, ellipse_y, ellipse_z, 'r-', LineWidth=2);

theta2 = linspace(0, 2*pi, 100);
circle_x = ((abs(cos(theta2)/R)).^2 + (abs(sin(theta1)/R)).^2).^(-1/2) .* cos(theta1) + 0.01;
circle_z = ((abs(cos(theta1)/R)).^2 + (abs(sin(theta1)/R)).^2).^(-1/2) .* sin(theta1) + R - 0.03;
circle_y = zeros(length(theta2)); % All points on z=0 plane
plot3(circle_x, circle_y, circle_z, 'b--', LineWidth = 2);

xlim([-0.2 0.2])
zlim([-0.05 0.1])
hold off

%% Assign pressure distribution

% Assign pressure to squares
pressure = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        if inside_superellipse(i, j)
            Pi(i, j) = CP_pressure.predictFcn([40, rad2deg(asin(X_rot(i,j)/R))]);
            y1(i,j) = b*(1-(abs(X_rot(i,j)/a).^n )^1/n);
            a1(i,j) = (1 - abs(Y_rot(i,j)/(y1(i,j))).^3).^2;
            P(i,j) = Pi(i, j) * a1(i,j);
        end
    end
end

%% Evaluate contact stress

% Define the target force and initial guess for Pm.
F_target = 2500;
Pm_initial = 125000;

% Define function 
fun = @(Pm) force(P, square_size, X, inside_superellipse, Pm, F_target);

% Use fzero to find the value of Pm.
P_mopt = fzero(fun, Pm_initial);

% Display the result.
fprintf('The value of Pm that produces a force of %.4f is %.4f\n', F_target, P_mopt);

function F = force(P, square_size, X, inside_superellipse, Pm, F_target)
    % Initialize the computed force.
    computedForce = 0;
    
    % Precompute maximum value of P for normalization.
    pmax = max(P(:));
    
    % Loop over the indices of X.
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            % Use the logical variable directly.
            if inside_superellipse(i, j)
                computedForce = computedForce + (P(i, j) / pmax) * Pm * square_size^2;
            end
        end
    end
    
    % Return the difference between the computed force and the target force.
    F = computedForce - F_target;
end

%% Plot contact stress

C = P/max(P(:))*P_mopt/1000;
figure;
surf(vertices_x, vertices_y, vertices_z, C, 'EdgeColor', 'none');

% Use colormap to indicate depth
colormap("parula");  % Jet colormap gives a good depth representation
cb = colorbar;  % Add color scale
ylabel(cb, "Normal stress (kPa)", 'Rotation',270)
caxis([0, max(C(:))]); % Scale color range to z-values

% Formatting the plot
shading interp; % Smooth color shading
grid on;
view(3); % Set 3D view
xlabel('X co-ordinate (m)');
ylabel('Y co-ordinate (m)');
zlabel('Z co-ordinate (m)');
axis equal
zlim([-0.08 0.0])
view([-45 45])
hold on
plot3(ellipse_x, ellipse_y, ellipse_z, 'r-', LineWidth=2);
hold off

%% Evaluate subsoil stess

% Define square area below the superellipse at depth
x_below = -3*a:square_size:3*a; % Square area grid (x-coordinates)
y_below = -3*a:square_size:3*a; % Square area grid (y-coordinates)
[X_below, Y_below] = meshgrid(x_below, y_below);

% Initialize vertical stress
stress = zeros(size(X_below));

di = depth*ones(size(vertices_z))+vertices_z;

% Calculate stress below the superellipse
for i = 1:size(X_below, 1)
    for j = 1:size(X_below, 2)
        sigma_z = 0; % Accumulated stress
        for m = 1:size(X, 1)
            for n = 1:size(X, 2)
                if inside_superellipse(m, n) && C(m, n) > 0
                    % Compute distance from evaluation point to surface square
                    r = sqrt((X_below(i, j) - X(m, n))^2 + (Y_below(i, j) - Y(m, n))^2);
                    % Contribution of each square
                    delta_sigma = square_size^2 * (3 * C(m, n) / (2 * pi * di(m, n)^2)) * 1 / ((r / di(m,n))^2 + 1)^5/2;
                    sigma_z = sigma_z + delta_sigma; % Accumulate stress
                end
            end
        end
        stress(i, j) = sigma_z; % Store result for this evaluation point
    end
end

%% Plot subsoil stess

% Plot vertical stress distribution below superellipse
figure;
surf(X_below, Y_below, stress, 'EdgeColor', 'none');
cb1 = colorbar;  % Add color scale
ylabel(cb1, "Normal stress (kPa)", 'Rotation',270)

xlabel('X co-ordinate (m)');
ylabel('Y co-ordinate(m)');
axis equal
xlim([-0.15 0.15])
ylim([-0.15 0.15])
view(2);
