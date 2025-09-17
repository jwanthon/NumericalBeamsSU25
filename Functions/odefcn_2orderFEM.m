%% odefcn_2orderFEM.m
%  Joseph Anthony
%
% Created:         5/15/25
% Last Modified:   9/17/25
%
% Description: Reduces the order of a system of ODEs in the form M.y'' = K.y
%   to a first order system. The first n functions define the actual n
%   solutions, while the second n functions are dummy functions
%
% INPUTS:
%   t: dummy variable used for the output
%   y: dummy variable used for the output
%   K: system's stiffness matrix
%   M: system's mass matrix
%   beta: linear damping (if unspecified, equals zero)
% OUTPUTS:
%   dydt: input to an ODE solver (i.e. ode45)

function dydt = odefcn_2orderFEM(t, y, K, M, beta)
    if ~exist('beta', 'var')
        beta = 0;
    end
    n = size(K, 1);
    invM = inv(M);
    dydt = zeros(2*n, 1);
    dydt(1:n) = y(n+1 : 2*n); % y'
    dydt(n+1:2*n) = invM*(K*y(1:n)-beta*dydt(1:n)); % y''
end

%% Explanation
% We know that M.y'' = K.y + By'. We need to decompose this into a system of
%   first-order ODEs. We will do this by implementing a solution vector, Y,
%   where Y(1:n) corresponds to y', and Y(n+1:2n) corresponds to y'.
%   Essentially, we will always be solving for the displacement and first
%   derivative, even though we ultimately only care about the displacement.
% 
% We also need a vector dydt which is 2n elements long: the first half is for the
%   derivatives of y, and the second half for the derivatives of y'.
%
% dydt(1:n) are the derivatives of y. So, dydy(1:n) = Y(n+1:2n).
% dydt(n+1:2n) are the derivatives of y'. We know that this is equal to
%   invM*(K.y+By'). y is the same as Y(1:n), while y' is the same as
%   dydt(1:n). We can multiply in our damping directly.