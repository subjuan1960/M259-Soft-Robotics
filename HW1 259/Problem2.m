N = 3:2:21;
dt = [1e0 1e-1 1e-2 1e-3];
%% simulation with N = 21 and dt = 0.01
N21 = N(end);
dt21 = dt(3);
[qlist, ulist, all_mid_v, all_mid_y, Nsteps] = p2sim(N21,dt21);

%% Plotting beam deformation in time march
figure(1);
for i = 1:Nsteps
    fprintf('Time = %f\n', (i-1) * dt21 );
    q = qlist(i,:);
    x = q(1:2:end);
    y = q(2:2:end);
    plot( x, y, 'ro-');
    title('Beam deformation');
    axis equal
    drawnow
end
%% Plotting mid-node position vs. Time
figure(2);
timeArray = (1:Nsteps) * dt21;
plot(timeArray, all_mid_y, 'k-');
title('Mid-node position vs. Time');
xlabel('Time, t [sec]');
ylabel('Position of mid-node, y [meter]');
%% Plotting mid-node velocity vs. Time
figure(3);
timeArray = (1:Nsteps) * dt21;
plot(timeArray, all_mid_v, 'k-');
title('Mid-node velocity vs. Time');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
%% Final deformed shape of the beam
qfinal = qlist(end,:);
xfinal = qfinal(1:2:end);
yfinal = qfinal(2:2:end);
figure(4)
plot( xfinal, yfinal, 'ro-');
title('Final Beam Deformation');
grid on
axis equal
%% Terminal velocity vs Number of Nodes
VterminalN = [];
for n = 3:2:21
    [~, ~, all_mid_vN, ~, ~] = p2sim(n,dt21);
    VterminalN = [VterminalN all_mid_vN(end)];
end

figure(5)
plot( N, VterminalN, 'ro-');
title('Terminal velocity vs Number of Nodes');
xlabel('Number of Nodes, N');
ylabel('Terminal Velocity of Nodes, v [meter/sec]');
grid on
%% Terminal velocity vs Varied timestep
Vterminaldt = [];
for i = 1:4
    [~, ~, all_mid_vdt, ~, ~] = p2sim(N21,dt(i));
    Vterminaldt = [Vterminaldt all_mid_vdt(end)];
end

figure(6)
semilogx( dt, Vterminaldt, 'ro-');
title('Terminal velocity vs Varied timestep');
xlabel('Timestep, [s]');
ylabel('Terminal Velocity of Nodes, v [meter/sec]');
grid on

