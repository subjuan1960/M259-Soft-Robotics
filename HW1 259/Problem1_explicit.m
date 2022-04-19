clear
clc
close all
%% Physical Parameters
% Number of vertices
N = 3;
% Time step size
dt = 1e-5;  
% Rod Length
RodLength = 0.1;
% discrete length
deltaL = RodLength / (N-1);
% Radius of Sphere
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;
% Density
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;
% Rod radius
r0 = 0.001;
% Young's Modulus
Y = 1e9; 
% gravity
g = 9.8;
% Viscosity
visc = 1000;
% Total time
totalTime = 10;
% Utility quantities
ne = N-1; % number of edges
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
% Geometry
nodes = zeros(N,2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
end
% Mass matrix
M = zeros(2*N, 2*N);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;
% Viscous Damping Matrix
C = zeros(6,6);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;
% Gravity
W = zeros(2*N, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;
% Initial DOF vector
q0 = zeros(2*N, 1);
for c = 1:N
    q0 (2 * c - 1 ) = nodes(c,1);
    q0 (2 * c) = nodes(c,2);
end
% New position and velocity
q = q0;
u = (q - q0) / dt;
curvature0 = 0;
% Number of time steps
Nsteps = round( totalTime / dt);
all_mid_y = zeros( Nsteps, 1);
all_mid_v = zeros( Nsteps, 1);
all_mid_y(1) = q(4);
all_mid_v(1) = u(4);

% calculates new q based on old q and plots in real time
for c = 2: Nsteps
    fprintf('Time = %f\n', (c - 1) * dt );
    q_new = q_new_calc(dt,q,u,curvature0,deltaL,EA,EI,M,W,C) % q_new calculation
    u = (q_new - q) / dt; % update u with u_new
    q = q_new; % update q with q_new
    figure(1);
    plot( q(1:2:end), q(2:2:end),  'ro-');
    axis equal
    drawnow

    % Store
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
end

figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t[sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
