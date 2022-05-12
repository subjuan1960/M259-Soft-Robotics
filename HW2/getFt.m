function [Ft, Jt] = getFt(x, refTwist)

% Ft is a vector of size 4N-1
% Jt is a square matrix of the same size

global GJ voronoiRefLen

ndof = numel(x);
ne = (ndof+1)/4 - 1; % number of edges

Ft = zeros(ndof, 1);
Jt = zeros(ndof, ndof);
for c = 2:ne
    ind = 4*c-7 : 4*c+3;
    node0 = [x(4*c-7), x(4*c-6), x(4*c-5)]; 
    node1 = [x(4*c-3), x(4*c-2), x(4*c-1)]; 
    node2 = [x(4*c+1), x(4*c+2), x(4*c+3)];
    theta_e = x(4*c-4);
    theta_f = x(4*c);
    l_k = voronoiRefLen(c);    
    [dF, dJ] = gradEt_hessEt(node0, node1, node2, ...
        theta_e, theta_f, refTwist(c), l_k, GJ);
    Ft(ind) = Ft(ind) - dF;
    Jt(ind, ind) = Jt(ind, ind) - dJ;
end

end
