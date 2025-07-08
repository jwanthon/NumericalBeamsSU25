%% odefcn_2orderFDM.m
%  Joseph Anthony
%
% Created:         5/15/25
% Last Modified:   7/8/25
%
% Description: Reduces the order of a system of ODEs in the form y'' = A.y
%   to a first order system. The first n functions define the actual n
%   solutions, while the second n functions are dummy functions
%
% INPUTS:
%   t: dummy variable used for the output
%   y: dummy variable used for the output
%   matrix: square matrix defining the system (i.e. FDM matrix)
%   beta: linear damping (if unspecified, equals zero)
% OUTPUTS:
%   dydt: input to an ODE solver (i.e. ode45)

function dydt = odefcn_2orderFDM(t, y, matrix, beta)
    if ~exist('beta', 'var')
        beta = 0;
    end
    n = size(matrix, 1);
    dydt = zeros(2*n, 1);
    dydt(1:n) = y(n+1 : 2*n); % y'
    dydt(n+1:2*n) = matrix * y(1:n) - beta*dydt(1:n); % y''
end
