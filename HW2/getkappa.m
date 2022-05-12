function kappa = getkappa( x, m1, m2 )
% Returns the curvature values at each node given the DOF vector, m1, and
% m2
nv = (numel(x) + 1 ) / 4;
ne = nv  - 1 ;

% Compute kappa
kappa = zeros(nv, 2);
for c=2:ne
    node0 = x(4*(c-2)+1:4*(c-2)+3); % (c-1)-th node
    node1 = x(4*(c-1)+1:4*(c-1)+3); % c-th node
    node2 = x(4*(c-0)+1:4*(c-0)+3); % (c+1)-th node
    m1e = m1(c-1,:); % m1 vector on (c-1)-th edge
    m2e = m2(c-1,:); % m2 vector on (c-1)-th edge
    m1f = m1(c,:); % m1 vector on (c)-th edge
    m2f = m2(c,:); % m2 vector on (c)-th edge
    
    kappaL = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f );
    
    kappa(c,1) = kappaL(1);
    kappa(c,2) = kappaL(2);
end

end
