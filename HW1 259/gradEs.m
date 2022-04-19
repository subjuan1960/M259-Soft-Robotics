function F = gradEs(xk, yk, xkp1, ykp1, l_k, EA)
%
% This function returns the derivative of stretching energy E_k^s with 
% respect to x_{k-1}, y_{k-1}, x_k, and y_k.
%
F = zeros(4,1);

F(1) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * xkp1 + 0.2e1 * xk);
F(2) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (-0.2e1 * ykp1 + 0.2e1 * yk);
F(3) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * xkp1 - 0.2e1 * xk);
F(4) = -(0.1e1 - sqrt((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) / l_k) * ((xkp1 - xk) ^ 2 + (ykp1 - yk) ^ 2) ^ (-0.1e1 / 0.2e1) / l_k * (0.2e1 * ykp1 - 0.2e1 * yk);

F = 0.5 * EA * l_k * F;

end
