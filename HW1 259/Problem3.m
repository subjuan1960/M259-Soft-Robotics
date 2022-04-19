%% Physical parameters
% Number of vertices
N = 50;

% Time step size
dt = 1e-3; % second

% Rod length
l = 1; % meter
d = 0.75; % load applied at

% Load
P = 20000; % Newtons (specified values of 2,000 and 20,000)

% Discrete length
deltaL = l / (N-1);
l_vec = [0:deltaL:l];
l_vec = abs(l_vec - d);
loc = min(l_vec);
x = find(l_vec == loc);

% Density
rho = 2700;

% Rod radius
R0 = 0.013;
r0 = 0.011;

% Young's modulus
Y = 70e9; % Using Y instead of E to avoid ambiguity

% Gravity
g = 9.8; % m/s^2

% Total time
totalTime = 1; % seconds

% Utility quantities
A = (R0^2 - r0^2)*pi;
I = (R0^4 - r0^4)*pi/4;
ne = N - 1; % Number of edges
EI = Y * I;
EA = Y * A;

% Geometry
nodes = zeros(N, 2);
for i = 1:N
    nodes(i,1) = (i-1) * deltaL;
end

% Mass matrix
m = pi*(R0^2-r0^2)*l*rho/(N-1);
M = eye(2*N,2*N)*m; % identity matrix times m

% Pressure location
Px = zeros(2*N,1);
Px(2*x) = -P;

% Initial DOF vector
q0 = zeros(2*N,1);
for i = 1:N
    q0 ( 2*i - 1 ) = nodes(i,1); % x coordinate
    q0 ( 2*i ) = nodes(i,2); % y coordinate
end

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% Number of time steps
Nsteps = round( totalTime / dt );
ymax = zeros(Nsteps, 1)

% Tolerance
tol = EI / l * 1e-3;



% Time marching scheme
for i=2:Nsteps
    
    fprintf('Time = %f\n', (i-1) * dt );
    
    q = q0; % Guess
    
    % Newton Raphson
    err = 10 * tol;
    while err > 1e-5
        % Inertia
        f = M / dt * ( (q-q0) / dt - u );
        J = M / dt^2;
        
        % Elastic forces
        %
        % Linear spring 1 between nodes 1 and N-1
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
            J(2*k-1:2*k+2,2*k-1:2*k+2) = ...
                J(2*k-1:2*k+2,2*k-1:2*k+2) + dJ;
        end
        
        % Bending spring between nodes 
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
                curvature0, deltaL, EI);
            f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
            J(2*k-3:2*k+2,2*k-3:2*k+2) = ...
                J(2*k-3:2*k+2,2*k-3:2*k+2) + dJ;
        end
        
        % Load
        f = f - Px;
        
        % Update
        q(3:2*N-1) = q(3:2*N-1) - J( ( 3:2*N-1 ),( 3:2*N-1 ) ) \ f( 3:2*N-1 );
        
        err = sum( abs( f( 3:2*N-1 )) );
    end

    % Update
    u = (q - q0) / dt; % Velocity
    q0 = q; % Old position
    
    % simulation of beam
    figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis([0 1 -0.3 0])
    drawnow
    grid on
    
    % Store
    ymax(i) = min(q);
end
%% Plotting max vertical displacement vs time
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, ymax, 'k-');
xlabel('Time, t [sec]');
ylabel('Maximum Vertical Displacement, y_m_a_x [meter]');
title('Maximum Vertical Displacement vs. Time');
grid on
%% Euler beam theory
c = min(d, l-d);
y_euler = (P*c*(l^2 - c^2)^1.5) / (9*sqrt(3)*EI*l);
fprintf('y_max from Euler beam theory is %1.4f m' ,y_euler);
%% Difference between simulation and Euler beam theory
max_disp = max(abs(ymax));
diff = abs(max_disp - y_euler);
