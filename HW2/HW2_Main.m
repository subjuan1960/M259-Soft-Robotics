%%
clear all
clc
% Global variables
global m dt unconsInd tol ScaleSolver maximum_iter
global Fg
global EI GJ voronoiRefLen kappaBar EA refLen

%% Parameters
nv = 100; % number of vertices
dt = 0.01; % Time step
RodLength = 0.2; % rod length
rho = 1000; % density
r0 = 1e-3; % cross sectional radius of the rod
Y = 10e6; % Young's modulus
nu = 0.5; % Poisson ratio
G = Y / (2*(1+nu)); % Shear modulus
natR = 0.02; % Natural radius of curvature
g = [0;0;-9.81];
tol = 1e-3; % Tolerance on force function (to be ...
% multiplied by characteristic bending force)
maximum_iter = 100; % maximum iterations in Newton's solver
totalTime = 5; % total simulation time

% Stiffness properties
EI = Y * pi * r0^4 / 4; % bending stiffness
EA = Y * pi * r0^2; % stretching stiffness
GJ = G * pi * r0^4 / 2; % twisting stiffness

% few other parameters
ndof = 4*nv - 1; % number of degrees of freedom
ne =  nv - 1; % number of edges
dm = (pi * r0^2 * RodLength) * rho / ne;
ScaleSolver = EI / RodLength^2; % Characteristic bending force

%% Geometry of the rod
nodes = zeros(nv, 3);
dTheta = RodLength / ne / natR;
for c=1:nv
    nodes(c,1) = natR * cos( (c-1) * dTheta );
    nodes(c,2) = natR * sin( (c-1) * dTheta );
end

%% Mass
m = zeros(ndof, 1);
for c=1:nv
    if c==1 || c==nv
        m( 4*(c-1) + 1: 4*(c-1) + 3 ) = dm/2;
    else
        m( 4*(c-1) + 1: 4*(c-1) + 3 ) = dm;        
    end        
end
for c=1:ne
    m( 4*c ) = dm/2 * r0^2; % I = 1/2*m*r^2
end

%% Gravity
garr = zeros(ndof, 1);
for c=1:nv
    garr( 4*(c-1) + 1: 4*(c-1) + 3 ) = g;
end
Fg = m .* garr; % Weight vector

%% Reference length
refLen = zeros(ne, 1);
for c=1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    refLen(c) = norm(dx);
end

%% Voronoi length
voronoiRefLen = zeros(nv, 1);
for c=1:nv
    if c==1
        voronoiRefLen(c) = 0.5*refLen(c);
    elseif c==nv
        voronoiRefLen(c) = 0.5*refLen(c-1);
    else
        voronoiRefLen(c) = 0.5*(refLen(c-1) + ...
            refLen(c));
    end
end

%% Reference directors
d1 = zeros(ne, 3); % Reference director, u
d2 = zeros(ne, 3); % Reference director, v

tangent = zeros(ne, 3);
for c=1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    tangent(c,:) = dx / norm(dx);
end

% Figure out a good choice for d1(1,:)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0, t1);
if abs(d1Tmp) < 1e-6
    t1 = [ 0 1 0 ];
    d1Tmp = cross(t0, t1);
end
d1(1,:) = d1Tmp / norm(d1Tmp);
d2Tmp = cross(t0, d1(1,:));
d2(1,:) = d2Tmp / norm(d2Tmp);

for c=2:ne
    t0 = tangent(c-1,:);
    t1 = tangent(c,:);
    d1_old = d1(c-1,:);
    d1_new = parallel_transport(d1_old, t0, t1);
    d1_new = d1_new / norm(d1_new);
    d1(c,:) = d1_new;
    d2_new = cross(t1, d1_new);
    d2(c,:) = d2_new / norm(d2_new);
end

%% Old and new degrees of freedom vector
x0 = zeros(ndof, 1); % Old dof (i.e. q0)
for c=1:nv
    x0( 4*(c-1) + 1 ) = nodes(c,1); % x-coord
    x0( 4*(c-1) + 2 ) = nodes(c,2); % y-cood
    x0( 4*(c-1) + 3 ) = nodes(c,3); % z-coord
end
x0(4:4:end) = 0; % theta

x = x0; % New dof (i.e. q)
u = (x - x0) / dt; % Old velocity vector

%% Constrained (fixed) and unconstrained (free) indices
consInd = 1:7; % Clamped; fixed_index
unconsInd = 8:ndof; % free_index

%% Material director
theta = x(4:4:end);
ne = numel(theta); % Number of edges
m1 = zeros(ne, 3);
m2 = zeros(ne, 3);

for c=1:ne
    cs = cos(theta(c));
    ss = sin(theta(c));
    d1_l = d1(c,:);
    d2_l = d2(c,:);
    m1_l = cs * d1_l + ss * d2_l;
    m1(c,:) = m1_l / norm(m1_l);
    m2_l = - ss * d1_l + cs * d2_l;
    m2(c,:) = m2_l / norm(m2_l);
end

%% Reference twist
refTwist = zeros(nv, 1);
refTwist = getRefTwist( d1, tangent, refTwist );

%% Natural curvature
kappaBar = getkappa( x, m1, m2 );

%% Simulation loop
Nsteps = round(totalTime/dt); % number of time steps

endZ = zeros(Nsteps, 1);

ctime = 0; % Current time

d1_old = d1;
d2_old = d2;
refTwist_old = refTwist;

for timeStep = 1:Nsteps % Main simulation step
    
    fprintf('t = %f\n', ctime);
    
    [x, d1, d2, refTwist] = objfun(x0, u, d1_old, d2_old, refTwist_old); % Main function that steps forward one time step

    % Plot
    x_coord = x(1:4:end);
    y_coord = x(2:4:end);
    z_coord = x(3:4:end);
    figure(1);
    clf();
    plot3(x_coord, y_coord, z_coord, 'ro-');
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    endZ(timeStep) = x(end);
    u = (x - x0) / dt; % Update velocity
    x0 = x; % New position becomes old position
    d1_old = d1; % New reference director becomes old reference director
    d2_old = d2;
    refTwist_old = refTwist; % New reference twist becomes old reference twist
    ctime = ctime + dt; % Current time
end

%%
figure(2);
time = (1:Nsteps) * dt;
plot( time, endZ, 'ro');
xlabel('Time, t [second]');
ylabel('Tip displacement, \delta_z [m]');
