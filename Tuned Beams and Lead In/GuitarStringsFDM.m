clc; clear all;
%% GuitarStringsFDM.m
%  Joseph Anthony
%
% Created:         5/5/25
% Last Modified:   5/7/25
%
% Description: Numerically solves the wave equation for a guitar string of
%  an abitrary length and determines the solution structure.
%
%% TODO:
% - Add audio conversion - maybe build as a function?
%
%% Building the FDM Matrix
% Let's model an undamped guitar string. 
% - Let u(x, t) be the string's displacement, with
%   utt being the second partial of u with respect to time, and uxx being the
%   second partial of u with respect to distance along the string. 
% - Let the string's tension be described by T [N] and the string's cross-sectional
%   density be described by rho [kg/m]. 
% - Finally, suppose the string is of a length L [m]. We must enforce
%   boundary conditions where u(0, t) and u(L, t) are both zero.

rho     = 1;                          % Parameter, density [kg/m^2]
T       = 1;                          % Parameter, tension [N]
L       = 1;                          % Parameter, length [m]

% Doing force analysis of the string, we see that:
%   rho * utt = T * uxx, or:
%   utt = T / rho * uxx                                          (eq. 1)
% 
% To figure out what our string should look like when it vibrates, we can
%   use the finite difference method (FDM) to determine the mode shapes of
%   the string.
% 
% To begin, let's create our mesh. That is, we subdivide the string's
%   length into equally sized chunks to perform analysis. We want there to
%   be n chunks.
% 
% Let the ith position on the string, xi, be equal to i * deltaX. i
%   will increment from 1 to n. We can describe this as a vector x. Since
%   the boundary conditions are zero, we only care about points internal to
%   the string.

n       = 100;                          % Parameter, mesh size

deltaX  = L / (n  + 1);
x       = linspace(deltaX, L - deltaX, n); 

% Using the wave equation, we know that utt depends only on a constant and
%   uxx. Knowing the difference quotient definition of the derivative, we
%   can approximate uxx using an algebraic system of equations of u(xi, t). 
% 
% The first derivative at xi, ux(xi, t), can be approximated by:
%   ux(xi, t) = (u(xi + deltaX, t) - u(xi - deltaX, t)) / deltaX
% 
% And by using another finite difference of these values at xi, we find:
%   uxx(xi, t) = ((uxi - deltaX, t) - 2(xi, t) + u(xi + deltaX, t)) /
%   deltaX^2                                                     (eq. 2)
% 
% But this is a system of equations depending only on values of xi! We can
%   thus describe this by the matrix-vector equation:
%   utt = -T / (rho * deltaX^2) * A.u                            (eq. 3)
%
%   where A is the finite difference matrix, and u is simply the displacement
%   across the length of the string, where the ith element describes the
%   displacement at xi. Note that the position information is baked into the
%   product A.u, so now our system is dependent only on time.
%
% Each row of the matrix should contain a series of 1, -2, 1 to describe
%   the computations done in the finite difference, where the -2 sits on
%   the main diagonal. Since the boundary conditions on the string are zero,
%   any 1s that would appear on the matrix in the topmost and bottommost
%   rows are also zeroes and do not occur in the system. Note that eq. 3
%   negates this matrix to have a positive main diagonal, so our matrix
%   A will be populated by -1s and 2s.
%
% Note that other FDM techniques (especially for ODEs) will place these
%   boundary conditions on the topmost and bottommost rows.

% Build the FDM matrix
temp1   = repmat(2, n, 1);
temp2   = repmat(-1, n-1, 1);
A       = diag(temp1) + diag(temp2, 1) + diag(temp2, -1);

clear temp1 temp2;

%% Mode Shapes
% Inspecting eq. 3, we can see that utt is defined by the matrix-vector
%   product A.u. Because of the way that this PDE operates, we can break
%   our solution into two pieces of information: the solution modes, and
%   the solution mode strengths.
% 
% Our boundary conditions inform the way that the matrix for uxx is built.
%   As a result, only certain solution shapes, or modes can be built on the
%   string. However, we don't know what specific modes will propogate until
%   our initial conditions in time (i.e., plucking the string) are
%   satisfied.
% 
% From eq. 3, we see that any solution vector u will be a scalar mutliple
%   of utt. That is, all solution vectors of A will be eigenvectors of the
%   system! We don't know how these eigenvectors will be scaled in the
%   solution without our time-domain initial conditions, but these would
%   describe the modes, or shapes.
% Note that eigenvectors oscilate in sign from element to element, and must
%   be sorted to remove this behavior.

modeCount          = 5;               % Parameter, # of displayed modes
[eigVecs, eigVals] = eig(A);

% Sort the eigenvectors and eigenvalues from smallest to largest mode
eigVecs = flip(eigVecs, 1);
eigVals = flip(flip(eigVals, 1),2);

% Remove oscillation from each eigenvector and plot each mode, from largest
%   to smallest strength
tiledlayout(1,2); 
nexttile; 
hold on;
for i = 1:modeCount
    temp1 = eigVecs(:,i);
    plot(temp1,"Marker",".", "LineStyle","none")
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position (i)');
ylabel('Relative Mode Strength');

%% Time-Domain Solutions
% Let's install a basic initial condition: that there is an initial
%   displacement of the string at some location xi, and that the string's
%   displacement will linearly interpolate between this location and each
%   endpoint.
% We can assume that the string will be plucked somewhere close to the
%   center of the string, as well. Significant 

initDisp = [1, 0.1 * L];             % Parameter, [strength, x location]

icent = round(initDisp(2)/deltaX);
initLine = [linspace(0,initDisp(1), icent-1),initDisp(1),linspace(initDisp(1),0, n-icent)];

% By inspecting our modes, we see that our solutions will be sinusoids that
%   with zeroes (or nodes) at the endpoints of the strings, and with
%   the ith mode having i+1 modes equidistant across the string.
% Using eq. 3, we can find that the angular frequency of each mode is of
%   the form:
%   wi = sqrt(T * lambdai / (rho * deltaX^2))           (eq. 3)
%   
%   where wi is the ith angular frequency [rad], and lambdai is the
%   eigenvalue associated with the ith eigenvector (or mode).
% Let's add this information to our graph.

% Determine angular frequencies
w = diag(sqrt(T*abs(eigVals)/(rho*deltaX^2)));

% Add frequency information to graph
for i = 1:modeCount
    legends(i) = string(sprintf('%.2f [rad/s]', w(n-i+1)));
end
legend(legends);
hold off;

% To solve our time domain solution, we know that the superposition of each
%   eigenvector, of difference strengths, will result in the initial
%   displacement. Let our eigenbasis be K, strengths c, and initial
%   condition u. Then K.c = u.

eigWeights = linsolve(eigVecs, initLine');

% The solution in the time domain, to these initial conditions, should
%   experience exponential decay. However, we can see what the frequency of
%   this string profile should sound like, with the given mode count.

solution = zeros(100,1);
for i = 1:modeCount
    solution = solution+eigVecs(:,i).*eigWeights(i);
end

nexttile
plot(solution)
title(sprintf('Solution with %d String Modes', modeCount));
xlabel('String Position (i)');
ylabel('Relative Strength');