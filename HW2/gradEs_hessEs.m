%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.

%%
function [dF, dJ] = ...
    gradEs_hessEs(node0, node1, ...
    l_k, EA)
%
% Inputs:
% node0: 1x3 vector - position of the first node
% node1: 1x3 vector - position of the last node
%
% l_k: reference length (undeformed) of the edge
% EA: scalar - stretching stiffness - Young's modulus times area
%
% Outputs:
% dF: 6x1  vector - gradient of the stretching energy between node0 and node 1.
% dJ: 6x6 vector - hessian of the stretching energy between node0 and node 1.

%% Gradient of Es
edge = (node1 - node0)'; % 3x1 edge vector
edgeLen = norm(edge);
tangent = edge / edgeLen;
epsX = edgeLen/l_k - 1;
dF_unit = EA * tangent * epsX;

dF = zeros(6,1);
dF(1:3) = - dF_unit;
dF(4:6) =   dF_unit;

%% Hessian of Es
Id3 = eye(3);
M = EA * ( ...
    (1/l_k - 1/edgeLen) * Id3 + ...
    1/edgeLen * (edge*edge')/ edgeLen^2 ...
    ); %Note edge * edge' must be 3x3

dJ = zeros(6,6);
dJ(1:3, 1:3) = M;
dJ(4:6, 4:6) = M;
dJ(1:3, 4:6) = - M;
dJ(4:6, 1:3) = - M;

end
