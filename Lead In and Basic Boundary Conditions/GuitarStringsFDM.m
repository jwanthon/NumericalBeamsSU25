clc; clear all;
addpath('Functions\');
%% GuitarStringsFDM.m
%  Joseph Anthony
%
% Created:         5/5/25
% Last Modified:   5/14/25
%
% Description: Numerically solves the wave equation for a guitar string of
%  an abitrary length using FDM and determines the solution structure.
%
% Required: 
%   solutionToAudio.m                                       (Function)
%   solutionToFFT.m                                         (Function)
%% TODO:
% - Add damped response
% - Frequency response
%   - Add audio conversion
%       - Maybe build as a function?
%   - Add graph of frequency response
%
%% Building the FDM Matrix
% Let's model an undamped guitar string. 
% - Let u(x, t) be the string's displacement, with
%   utt being the second partial of u with respect to time, and uxx being the
%   second partial of u with respect to distance along the string. 
% - Let the string's tension be described by T [N] and the string's cross-sectional
%   density be described by ρ [kg/m]. 
% - Finally, suppose the string is of a length L [m]. We must enforce
%   boundary conditions where u(0, t) and u(L, t) are both zero.
%
% Doing force analysis of the string, we see that:
%   ρ * utt = T * uxx, or:
%   utt = T / ρ * uxx                                          (eq. 1)

rho = 1;                          % Parameter, density (ρ) [kg/m^2]
T   = 1;                          % Parameter, tension (T) [N]
L   = 1;                          % Parameter, length  (L) [m]

% We will numerically solve for the string's vibration using the finite difference
%   method (FDM). To begin, let's create a mesh and subdivide the string's
%   length into n equally sized chunks.

n   = 100;                        % Parameter, mesh size

% We want the ith position of the string, xi, to be equal to i * ∆x. We also
%   only care about points internal to the string, since we know by the
%   boundary conditions that the points on the ends of the strings are
%   stationary!

deltax  = L / (n  + 1);

% Using FDM on the wave equation allows us to create an algebraic system of
%   equations based on approximating uxx as a finite differnence of u(xi,
%   t) for fixed t.
% 
% The first derivative at xi, ux(xi, t), can be approximated by:
%   ux(xi) = (u(xi+∆x, t) - u(xi-∆x, t)) / ∆x
% 
% And by using another finite difference of these values at xi, we find:
%   uxx(xi, t) = (u(xi-∆x, t) - 2u(xi, t) + u(xi+∆x, t)) / ∆x^2 (eq. 2)
%
% Hence, uxx depends only on the values of u at different points on the
%   mesh. That is, uxx(xi) = u(x(i-1)) - 2u(xi) + u(x(i+1)).
%   Our boundary conditions enforce special conditions for i = 1 and i = n.
%   When i = 1, the point to the left of xi is an end of the string, where
%   u = 0, and:
%   uxx(x1) = (0 - 2u(x1, t) + u(x2, t)) / ∆x^2 
%
% Similarly, when i = n, the point to the right of xi is another
%   end of the string, where u = 0, and:
%   uxx(xn) = (u(x(n-1), t) - 2u(xn, t) + 0) / ∆x^2
% 
% That is, uxx can be represented by a system of linear equations, all 
%   depending only on points internal to the string. This can also be
%   expressed by the matrix-vector equation:
%
%   utt = -T / (ρ * ∆x^2) * A.u                            (eq. 3)
%
%   where A is the finite difference matrix, and u is simply the
%   displacement across the length of the string for all values xi at time
%   t. Our matrix is in the following form:
%
%    2      -1        0     ...  0
%   -1       2       -1     ...  0
%    0      -1        2     ...  0
%    .       .        .   .      .
%    .       .        .      .  -1
%    0       0        0   -1     2
%
%   with 2s on the main diagonal, and -1s directly above and below the main
%   diagonal.

% Build the FDM matrix
temp1   = repmat(2, n, 1);
temp2   = repmat(-1, n-1, 1);
FDM       = diag(temp1) + diag(temp2, 1) + diag(temp2, -1);

clear temp1 temp2 % Cleanup
disp('FDM Matrix Built');
%% Mode Shapes
% From eq. 3, we know that utt is defined by both the position information,
%   which now depends only on time, and the FDM matrix A. Note that A is
%   built using only our boundary conditions! 
% 
% We can see that any basic solution shape, or mode, of the system, when
%   passed through A, will output a scalar multiple of itself. That is, the
%   modes of the system on the string are the eigenvectors of A.

modeCount          = 5;               % Parameter, # of displayed modes
[eigVecs, eigVals] = eig(FDM);

% Sort the eigenvectors and associated eigenvalues from smallest to largest
eigVecs = flip(eigVecs, 1);
eigVals = flip(flip(eigVals, 1),2);

% Plot the specified number of modes
fig2 = tiledlayout(1,2); 
nexttile; 
hold on;
for i = 1:modeCount
    plot(eigVecs(:,i),"Marker",".", "LineStyle","none");
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position (i)');
ylabel('Displacement [m]');

% By inspecting our modes, we see that our solutions will be sinusoids that
%   with zeroes (or nodes) at the endpoints of the strings, and with
%   the ith mode having i-1 nodes equidistant across the string.
% Using eq. 3, we can find that the angular frequency of each mode is of
%   the form:
%   wi = sqrt(T * λi / (ρ * ∆x^2))           (eq. 4)
%   
%   where wi is the ith angular frequency [rad], and λi is the eigenvalue 
%   associated with the ith eigenvector (or mode).

% Determine angular frequencies
angFreq = diag(sqrt(T*abs(eigVals)/(rho*deltax^2)));

% Add frequency information to graph
for i = 1:modeCount
    legends(i) = string(sprintf('%.2f [rad/s]', angFreq(n-i+1)));
end
legend(legends);
hold off;
clear legends i % Cleanup

%% Time-Domain Solutions
% Let's install some initial conditions on the string at time t = 0.

% Create initial conditions at one point
pointPos = [5, 0.5 * L];             % Parameter, [displacement [m], location]
pointVel = [-4, 0.5 * L];            % Parameter, [velocity [m/s], location]

% Discritize initial conditions, and assume position and velocity changes
%   linearly across the string.
centerpoint = round(pointPos(2)/deltax);
initPos = [linspace(0,pointPos(1), centerpoint-1),pointPos(1),linspace(pointPos(1),0, n-centerpoint)];
centerpoint = round(pointVel(2)/deltax);
initVel = [linspace(0,pointVel(1), centerpoint-1),pointVel(1),linspace(pointPos(1),0, n-centerpoint)];

clear pointPos pointVel % Cleanup

% To solve our time domain solution, we need to solve eq. 3. From ODEs, we 
%   know that all of the eigenvectors will be sines and cosines in the time 
%   domain, so we will know that for the ith element of the solution vector,
%   ui = ki(ai * cos(wi * t) + bi * sin(wi * t)).
%   
%   We can solve for the vectors a and b using our initial
%   conditions. At t = 0, ui = ai*ki, and ui' = ki*wi*bi.

scaledFDM = -T / (rho * deltax^2) * FDM;      % Matrix in eq. 3
aCoeff = linsolve(scaledFDM, initPos');       % Solving K.a = u
bCoeff = linsolve(scaledFDM, initVel');       % Solving k.(b.*w) = u'
bCoeff = bCoeff./angFreq;                     % Dividing out w from b.*w

% Solve for the low-resolution solution (animated)
i = 1;
timemax = 10;
timestep = 0.02;
for t = 0:timestep:timemax
    for j = 1:modeCount
        anim_solution(:,i) = aCoeff.*eigVecs(:,j).*cos(t*angFreq) + bCoeff.*eigVecs(:,j).*sin(t*angFreq);
    end
    i = i+1;
end
length = i -1;

% Solve for high-res solution (used in FFT)
i = 1;
for t = 0:1/10000:1
    for j = 1:modeCount
        solution(:,i) = aCoeff.*eigVecs(:,j).*cos(t*angFreq) + bCoeff.*eigVecs(:,j).*sin(t*angFreq);
    end
    i = i+1;
end

%% Animated Solution
% Plot an animated graph of the solution
[solution_min, solution_max] = bounds(anim_solution);    % Determine axes
solution_min = min(solution_min);
solution_max = max(solution_max);

nexttile;

while 1
    for i = 1:length
    plot(anim_solution(:,i))
    title(sprintf('Undamped Solution with %d String Modes', modeCount));
    xlabel('String Position (i)');
    ylabel('Relative Strength');
    axis([1 n solution_min solution_max]);
    drawnow
    pause(timestep)
    end
end
%% 
[transform, freqspace] = solutionToFFT(solution, 50, 1/10000);
plot(freqspace, transform);
%% References
% [1] A. Struthers. (2025). Mathematical modeling of a whammy bar kalimba
%   [Mathematica slides]. Available:
%   https://github.com/AllanStruthersMTU/Spring-Classes-2025/blob/main/ComplexAnalysis/April25LakeStateKalimba.nb.