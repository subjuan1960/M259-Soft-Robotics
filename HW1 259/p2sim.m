function [qlist, ulist, all_mid_v, all_mid_y, Nsteps] = p2sim(N,dt)

    % Rod length
    RodLength = 0.1; % meter
    
    % Discrete length
    deltaL = RodLength / (N-1);
    
    % Radius of spheres
    R = deltaL/10;
    R_mid = 0.025;
    
    
    % Density
    rho_metal = 7000;
    rho_f = 1000;
    rho = rho_metal - rho_f;
    
    % Rod radius
    r0 = 0.001;
    
    % Young's modulus
    Y = 1e9; % Using Y instead of E to avoid ambiguity
    
    % Gravity
    g = 9.8; % m/s^2
    
    % Viscosity
    visc = 1000; % Pa-s
    
    % Total time
    totalTime = 50; % seconds
    
    % Utility quantities
    ne = N - 1; % Number of edges
    EI = Y * pi * r0^4 / 4;
    EA = Y * pi * r0^2;
    
    % Geometry
    nodes = zeros(N, 2);
    for i = 1:N
        nodes(i,1) = (i-1) * deltaL;
    %     nodes(c,2) = 0;
    end
    
    % Mass matrix
    M = zeros(2*N,2*N);
    for i = 1:2N
        if i == N || i == N+1
        M(i,i) = 4/3*pi*R_mid^3*rho_metal;
        else
        M(i,i) = 4/3*pi*R^3*rho_metal;
        end
    end
    
    % Viscous damping matrix
    C = zeros(2*N,2*N);
    for i = 1:2*N
        if i == N || i == N+1
        C(i,i) = 6*pi*visc*R_mid;
        else
        C(i,i) = 6*pi*visc*R;
        end
    end
    
    
    % Gravity
    W = zeros(2*N,1);
    for i = 1:N
        if i == (N+1)/2
            W(2*i) = -4/3*pi*R_mid^3*rho*g;
        else
            W(2*i) = -4/3*pi*R^3*rho*g;
        end
    end
    
    % Initial DOF vector
    q0 = zeros(2*N,1);
    for i=1:N
        q0 ( 2*i - 1 ) = nodes(i,1); % x coordinate
        q0 ( 2*i ) = nodes(i,2); % y coordinate
    end
    
    % New position and velocity
    q = q0; % DOF vector
    u = (q - q0) / dt; % Velocity vector

    % lists for storing position and velocity
    qlist = [];
    ulist = [];
    qlist = q';
    ulist = u';
    
    % Number of time steps
    Nsteps = round( totalTime / dt );
    all_mid_y = zeros( Nsteps, 1); % y-position of R_mid
    all_mid_v = zeros( Nsteps, 1); % y-velocity of R_mid
    
    all_mid_y(1) = q(N+1);
    all_mid_v(1) = u(N+1);
    
    % Tolerance
    tol = EI / RodLength^2 * 1e-3;
    
    % Time marching scheme
    for i=2:Nsteps
        
        fprintf('Time = %f\n', (i-1) * dt );
        
        q = q0; % Guess
        
        % Newton Raphson
        err = 10 * tol;
        while err > tol
            % Inertia
            f = M / dt * ( (q-q0) / dt - u );
            J = M / dt^2;
            
            %
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
            
            % Viscous force
            f = f + C * ( q - q0 ) / dt;
            J = J + C / dt;
            
            % Weight
            f = f - W;
            
            % Update
            q = q - J \ f;
            
            err = sum( abs(f) );
        end
    
        % Update
        u = (q - q0) / dt; % Velocity
        q0 = q; % Old position
        
%         % Plotting beam deformation in time march
%         figure(1);
%         title('Beam deformation');
%         plot( q(1:2:end), q(2:2:end), 'ro-');
%         axis equal
%         drawnow
        
        
        % Store
        qlist = [qlist; q'];
        ulist = [ulist; u'];
        all_mid_y(i) = q(N+1);
        all_mid_v(i) = u(N+1);
    end
end