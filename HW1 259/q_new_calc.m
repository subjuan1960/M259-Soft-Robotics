function q_new = q_new_calc(dt,q,u,curvature0,deltaL,EA,EI,M,W,C)
    
    dEs1 = gradEs(q(1), q(2), q(3), q(4), deltaL, EA);
    dEs2 = gradEs(q(3), q(4), q(5), q(6), deltaL, EA);
    dEb = gradEb(q(1), q(2), q(3), q(4), q(5), q(6), ...
            curvature0, deltaL, EI);
    dEs1 = [dEs1; 0; 0];
    dEs2 = [0; 0; dEs2];

    dE = dEs1 + dEs2 + dEb;
    q_new = dt * (inv(M) * dt * (W - dE - C*u) + u) + q;
end