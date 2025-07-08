clc; clear all;
addpath('Functions\');
%% BIBIBeamFDM.m
%  Joseph Anthony
%
% Created:         5/14/25
% Last Modified:   5/15/25
%
% Description: Numerically solves the Euler-Bernoulli beam equation 
%   for built in-built in beam using FDM and determines the solution 
%   structure.
%
%% TODO:
% - Solve ODE correctly
%% Building the FDM Matrix
% We will now be working with the dynamic Euler-Bernoulli beam equation: 
%   ρ * utt = Γ * uxxx
%   which assumes the beam is homogenous. We say that Γ is the flexural
%   rigidity and ρ is the mass per unit length.

rho         = 10;        % Mass per unit length [kg/m]
gamma       = 1;        % Flexural rigidity [N·m²]
L           = 100;        % Beam length [m]
n           = 100;      % Mesh size
modeCount   = 5;        % Number of displayed modes

% We will be using a fourth-order central finite difference approximation 
%   for the RHS, where we are only concerned with points internal to the 
%   beam, for constant t:
%   uxxxx ≈ u(xi - 2Δx) - 4u(xi - Δx) + 6u(xi) - 4u(xi + Δx) + u(xi + 2Δx)
%
% For our BCs, we will say that our beam is built in (BI), or clamped at
%   both sides. This means that each side experiences no displacement, 
%   u(0, t) = u(L, t) = 0, and no moment, u'(0, t) = u'(L, t) = 0. We can
%   install this in our FDM matrix.
%
% We know that:
%   uxxxx(1)  = u(-1) - 4u(0) + 6u(1) - 4u(2) + u(3)
%   uxxxx(2)  =          u(0) - 4u(1) + 6u(2) - 4u(3) + u(4)
%   ux(0) = 0 = u(-1) + 2u(0) -  u(1)
% But u(0) = 0, so by eq. 3 u(-1) = u(1). Hence, the upper corner of our
%   matrix will be:
%
%   7   -4  1   0   0  ...  % The upper left corner of the matrix is a 7,
%   -4  6   -4  1   0  ...  % with 6's on the main diagonal with bands of
%   1   -4  6   -4  1  ...  % -4's and 1's radiating off of it. The bottom-
%   0   1   -4  6   -4 ...  % right corner has the same structure.

% Building the FDM matrix
FDM = diag(repmat(-4, n-1, 1), 1) + diag(ones(n-2, 1), 2);
FDM = FDM + FDM';
FDM = FDM + diag(repmat(6, n, 1));
FDM(1,1) = 7;
FDM(n,n) = 7;
disp('FDM Matrix Built');

%% Mode Shapes
% Now our system is in the form:
%   ρ * utt = Γ * A.u
%   where A is our FDM matrix. We can once again determine our solution
%   modes from the eigenvectors of A. Note that these mode shapes are not
%   sinusoids like for the guitar string!

% Find and sort eigenvectors and eigenvalues of the FDM matrix
[eigVecs, eigVals] = eig(FDM);
eigVecs = flip(eigVecs, 1);
eigVals = flip(flip(eigVals, 1),2);

% Plot the specified number of modes

fig1 = tiledlayout(1,2); 
nexttile; 
hold on;
for i = 1:modeCount
    plot(eigVecs(:,i),"Marker",".", "LineStyle","none");
end
title(sprintf('First %d Mode Shapes of a BI-BI Beam', modeCount));
xlabel('Beam Position (i)');
ylabel('Displacement [m]');
hold off

%% Time-Domain Solution
% Like with the guitar string, we now have a second-order system of ODEs.
%   This means that we can solve for its time-domain solution. Since the
%   modes are not sines and cosines, we can't solve it directly: we must
%   solve it numerically.

% ICs
% We will assume that our loading force will be at the centerpoint of the
%   beam, and that there is no initial velocity on the beam.

initDisp = 1;          % Initial displacement at beam center [m]
timeMax  = 100;        % End point of simulation
timestep = .01;        % Timestep

IC = zeros(2*n, 1);
IC(1:n) = BIBIBeamICs(n, L, gamma, initDisp);
scaledFDM = FDM * gamma / rho;

[t,y] = ode45(@(t,y) odefcn_2orderFDM(t, y, scaledFDM), 0:timestep:timeMax, IC);
nexttile;
surf(1:n,t,y(:,1:n), "LineStyle","none");
camlight('headlight');
title('Beam Displacement over Time');
xlabel('Beam Position (i)');
ylabel('Time [sec]');
zlabel('Displacement [m]');

fig2 = visual_sparseMatrix(FDM);
fig2.Title = "FDM Matrix";

%% Animated Solution
% % Plot an animated graph of the solution
% [solution_min, solution_max] = bounds(solution);    % Determine axes
% solution_min = min(solution_min);
% solution_max = max(solution_max);
% 
% nexttile;
% while 1
%     for i = 1:size(solution,1)
%         plot(solution(:,i));
%         title(sprintf('Undamped Solution with %d Beam Modes', modeCount));
%         xlabel('Beam Position (i)');
%         ylabel('Relative Strength');
%         axis([1 n solution_min solution_max]);
%         drawnow
%         pause(.01);
%     end
% end