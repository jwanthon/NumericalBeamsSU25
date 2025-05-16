%% BIBIBeamICs.m
%  Joseph Anthony
%
% Created:         5/13/25
% Last Modified:   5/13/25
%
% Description: Uses displacement and velocity at the centerpoint of a
%   fixed-fixed beam to determine the magnitude of deflection and velocity
%   across the beam.
%
% INPUTS:
%   n:      mesh size
%   length: beam length
%   gamma:  flexural rigidity 
%   disp:   initial displacement at centerpoint
%   vel:    initial velocity at centerpoint
% OUTPUTS:
%   dispV:  n-length displacement vector

function vec_disp = BIBIBeamICs(n, length, gamma, disp)
    % Backsolve for force
    force = disp * 192 * gamma / length^3;
    vec_disp = zeros(n,1);
    % Populate displacement vector
    dx = length/(n+1);
    for i = 1:n/2
        vec_disp(i) = -force/48/gamma*(i*dx)^2*(3*length-4*i*dx);
    end
    vec_disp = vec_disp + flip(vec_disp);

    % Correct duplicate in the middle if n is odd
    if mod(n, 2) == 1
            vec_disp(n/2) = 1/2 * vec_disp(n/2);
    end
end

%% References
% [1] "Beam Deflection Tables", MechaniCalc. Accessed May 14, 2025.
%   [Online]. Available:https://mechanicalc.com/reference/beam-deflection-tables