%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.

%%
function [dF, dJ] = ...
    gradEt_hessEt(node0, node1, node2, ...
    theta_e, theta_f, refTwist, ...
    l_k, GJ)
%
% Inputs:
% node0: 1x3 vector - position of the node prior to the "twisting" node
% node1: 1x3 vector - position of the "twisting" node
% node2: 1x3 vector - position of the node after the "twisting" node
%
% theta_e: scalar - twist angle of the first edge
% theta_f: scalar - twist angle of the second (last) edge
%
% l_k: voronoi length (undeformed) of the turning node
% GJ: scalar - twisting stiffness
%
% Outputs:
% dF: 11x1  vector - gradient of the twisting energy at node1.
% dJ: 11x11 vector - hessian of the twisting energy at node1.

%% Computation of gradient of the twist
gradTwist = zeros(11,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

gradTwist(1:3) = -0.5 / norm_e * kb;
gradTwist(9:11) = 0.5 / norm_f * kb;
gradTwist(5:7) = -(gradTwist(1:3)+gradTwist(9:11));
gradTwist(4) = -1;
gradTwist(8) = 1;

%% Computation of hessian of twist
DDtwist = zeros(11, 11);

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;

D2mDe2 = -0.25 / norm2_e * ( kb' * (te + tilde_t) ...
    + (te + tilde_t)' * kb);
D2mDf2 = -0.25 / norm2_f  * ( kb' * (tf + tilde_t) ...
    + (tf + tilde_t)' * kb );
D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - ...
    kb' * tilde_t );
D2mDfDe = D2mDeDf';

DDtwist(1:3,1:3) = D2mDe2;
DDtwist(1:3, 5:7) = -D2mDe2 + D2mDeDf;
DDtwist(5:7, 1:3) = -D2mDe2 + D2mDfDe;
DDtwist(5:7, 5:7) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
DDtwist(1:3, 9:11) = -D2mDeDf;
DDtwist(9:11, 1:3) = -D2mDfDe;
DDtwist(9:11, 5:7) = D2mDfDe - D2mDf2;
DDtwist(5:7,9:11) = D2mDeDf - D2mDf2;
DDtwist(9:11,9:11) = D2mDf2;
    
%% Gradient of Et
integratedTwist = theta_f - theta_e + refTwist;
dF = GJ/l_k * integratedTwist * gradTwist;

%% Hessian of Eb
dJ = GJ/l_k * (integratedTwist * DDtwist + gradTwist*gradTwist');

end
