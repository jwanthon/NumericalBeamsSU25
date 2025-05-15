%% odefcn_2order.m
%  Joseph Anthony
%
% Created:         5/15/25
% Last Modified:   5/15/25
%
% Description: Reduces the order of a system of ODEs in the form y'' = A.y
%   to a first order system. The first n functions define the actual n
%   solutions, while the second n functions are dummy functions
%
% INPUTS:
%   t: dummy variable used for the output
%   y: dummy variable used for the output
%   matrix: square matrix defining the system (i.e. FDM matrix)
% OUTPUTS:
%   dydt: input to an ODE solver (i.e. ode45)

function dydt = odefcn_2order(t, y, matrix)
    n = size(matrix, 1);
    dydt = zeros(2*n, 1);
    dydt(1:n) = y(n+1 : 2*n);
    dydt(n+1:2*n) = matrix * y(1:n);
end
