clear
clc
close all
%% physical parameters
% Number of vertices
N = 3;

% Time step size
dt = 0.01;

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

% Number of time steps
Nsteps = round( totalTime / dt);
all_mid_y = zeros( Nsteps, 1);
all_mid_v = zeros( Nsteps, 1);

all_mid_y(1) = q(4);
all_mid_v(1) = u(4);

% Tolerance
tol = EI / RodLength^2 * 1e-3;
% Time marching scheme
for c = 2 : Nsteps
    fprintf('Time = %f\n', (c - 1) * dt );

    q = q0; % guess

    % Newton Raphson
    err = 10 * tol;
    while err > tol
        % Inertia
        f = M / dt * ( (q - q0) / dt - u);
        J = M / dt^2;
        % Elastic forces
        % Linear spring 1 between nodes 1 and 2
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(1:4) = f(1:4) + dF;
        J(1:4, 1:4) = J(1:4, 1:4) + dJ;
    
        % Linear spring 2 between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(3:6) = f(3:6) + dF;
        J(3:6, 3:6) = J(3:6, 3:6) + dJ;
    
        % Bending spring between nodes 1, 2 and 3
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0;
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
            curvature0, deltaL, EI);
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
            curvature0, deltaL, EI);
        f(1:6) = f(1:6) + dF;
        J(1:6, 1:6) = J(1:6, 1:6) + dJ;
    
        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;
    
        % Weight
        f = f - W;
    
        % update
        q = q - J \ f;
    
        err = sum( abs(f));
    end

    % Calculate q using N-R

    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position

    figure(1);
    plot( q(1:2:end), q(2:2:end),  'ro-');
    axis equal
%     xlim([0 0.1]);
%     ylim([-0.1 0]);
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


