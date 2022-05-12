function [Fs, Js] = getFs(x)

% Fs is a vector of size 4N-1
% Js is a square matrix of the same size

global EA refLen

ndof = numel(x);
ne = (ndof+1)/4 - 1; % number of edges

Fs = zeros(ndof, 1);
Js = zeros(ndof, ndof);
for c = 1:ne
    ind = [4*c-3, 4*c-2, 4*c-1, 4*c+1, 4*c+2, 4*c+3];
    node0 = [x(4*c-3), x(4*c-2), x(4*c-1)]; 
    node1 = [x(4*c+1), x(4*c+2), x(4*c+3)]; 
    l_k = refLen(c);    
    [dF, dJ] = gradEs_hessEs(node0, node1, l_k, EA);
    Fs(ind) = Fs(ind) - dF;
    Js(ind, ind) = Js(ind, ind) - dJ;
end

end
