function [Fb, Jb] = getFb(x, m1, m2)

% Fb is a vector of size 4N-1
% Jb is a square matrix of the same size

global EI voronoiRefLen kappaBar

ndof = numel(x);
ne = (ndof+1)/4 - 1; % number of edges

Fb = zeros(ndof, 1);
Jb = zeros(ndof, ndof);
for c = 2:ne
    ind = 4*c-7 : 4*c+3;
    node0 = [x(4*c-7), x(4*c-6), x(4*c-5)]; 
    node1 = [x(4*c-3), x(4*c-2), x(4*c-1)]; 
    node2 = [x(4*c+1), x(4*c+2), x(4*c+3)]; 
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);
    l_k = voronoiRefLen(c);    
    [dF, dJ] = ...
        gradEb_hessEb(node0, node1, node2, ...
        m1e, m2e, m1f, m2f, kappaBar(c,:), l_k, EI);    
    Fb(ind) = Fb(ind) - dF;
    Jb(ind, ind) = Jb(ind, ind) - dJ;
end

end
