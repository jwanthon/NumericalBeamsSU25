clc; clear all;
addpath('Functions\');
%% NumericFEMTestbench.m
%  Joseph Anthony
%
% Created:          6/24/25
% Last Modified:    6/24/25
%
% Description: Testbench for different properties of and ways to generate 
%   FEM stiffness matrices using numeric methods. This includes strictly
%   numerical ways of generating basis functions, testing the conditioning
%   of different numerical stiffness matrices, and seeing the range of
%   feasibility for numerical basis functions.
%
%% Numerical FEM Generation
tic
clc; clear all;
n = 1000;         % Mesh size, number of interior points
shapes = 100;     % Number of shape functions, typically equal to n
L = 1;          % Beam length
modeCount = 5;  % Number of displayed modes

xvals       = linspace(0,L,n+2);
beamBasis   = zeros(shapes, n+2);
DbeamBasis  = zeros(shapes,n+1);
D2beamBasis = zeros(shapes,n);
deltax      = L/(n+1);

% Generate shape functions
%   Equally spaced nodes x(x-L)(x-L/2) ...
% for i = 1:shapes
%     beamBasis(i,:) = xvals.*xvals.*(xvals - L * ones(1, n+2)).^2;
%     for j = 1:i
%         beamBasis(i,:) = beamBasis(i,:).*(xvals-j*L/i*ones(1,n+2));
%     end
% end
%   Chebyshev Polynomials of the First Kind
% for i = 1:shapes
%     beamBasis(i,:) = chebyshevT(2*i-1, xvals).*(xvals-L*ones(1,n+2));
% end
%   Sin(ix)
for i = 1:shapes
    beamBasis(i,:) = xvals*sin(pi*i*xvals/L);
end
%   Tents
% for i = 1:shapes
%     cent = i*L/(shapes+1);
%     dist   = L/(shapes+1);          % Distance on the side of each tent function 
%     for j = 1:length(xvals)
%         if (xvals(j) >= cent - dist && xvals(j) <= cent)
%             beamBasis(i, j) = dist*(xvals(j)-cent+dist);
%         elseif (xvals(j) >= cent && xvals(j) <= dist+cent)
%             beamBasis(i,j) = -dist*(xvals(j)-cent-dist);
%         end
%     end
% end
% beamBasis(:, 1) = zeros(shapes, 1);
% beamBasis(:, end) = beamBasis(:, 1);

% Normalize and scale basis functions
for i = 1:shapes
    beamBasis(i,:) = beamBasis(i,:)/max(abs(beamBasis(i,:)));
    beamBasis(i,:) = beamBasis(i,:);
end

% See shape functions
% hold on
% for i = 1:shapes
%     plot(xvals,beamBasis(i,:));
% end
% return

% Generate second derivatives
for i = 1:shapes
    DbeamBasis(i,:) = diff(beamBasis(i,:))/deltax;
    D2beamBasis(i,:) = diff(DbeamBasis(i,:))/deltax;
end

% Build stiffness matrix
K = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = DbeamBasis(row,:).*DbeamBasis(col,:);
        K(row,col) = trapz(product);
    end
end
K = K + K' - diag(diag(K));

% Find and sort eigenbasis from lowest mode to highest
[eigVecs, eigVals] = eig(K);
[d, index] = sort(diag(abs(eigVals))); % Sort based on magnitude
eigVals = eigVals(index,index);
eigVecs = eigVecs(:, index);

% Display eigenvalues
% plot(log(diag(abs(eigVals))), 'Marker','o', 'LineStyle', 'none');
% title(sprintf('Stiffness Matrix Eigenvalues for %d Shapes', shapes))
% ylabel('Eigenvalue Magnitude (Log)');
% % loglog(diag(abs(eigVals)), 'Marker','o', 'LineStyle', 'none');
% return

% Create the mode shapes by scaling the basis functions by each eigenvector
modeShapes = beamBasis'*eigVecs;
modeShapes = modeShapes';
disp('Mode shapes calculated.');

% Normalize eigenfunctions
for i = 1:shapes
    modeShapes(i,:) = modeShapes(i,:)/max(abs(modeShapes(i,:)));
end

figure()
% Plot modecount eigenvectors
hold on;
for i = 1:modeCount
     plot(xvals,modeShapes(i,:),"Marker",".");
end
title(sprintf('First %d Mode Shapes of a Guitar String', modeCount));
xlabel('String Position [m]');
ylabel('Relative Displacement');
hold off

toc
% Graph Eigenvalues
figure()
eigplot = diag(eigVals);
eigplot = eigplot/eigplot(1);
rooteigplot = sqrt(eigplot');
plot(rooteigplot, "LineStyle", "none", "Marker", ".");
title("Square Roots of Eigenvalues");
ylabel("sqrt(λ)")
%% Conditioning Testing
% Finds the condition number of different sizes of FEM matrices using
%   different shape function generators.
clc; clear all;

n = 1000;       % Mesh size, number of interior points
shapes = 100;    % Number of shape functions, typically equal to n
L = 1;          % Beam length
modeCount = 5;  % Number of displayed modes

xvals       = linspace(0,L,n+2);
beamBasis   = zeros(shapes, n+2);
DbeamBasis  = zeros(shapes,n+1);
D2beamBasis = zeros(shapes,n);
deltax      = L/(n+1);
condnumbs_nodes = zeros(1, shapes);
condnumbs_cheby = zeros(1, shapes);


% Find conditioning for equally spaced zeros
K = zeros(shapes);

for i = 1:shapes

    % Equally spaced zeros
    beamBasis(i,:) = xvals.*xvals.*(xvals - L * ones(1, n+2));
        for j = 1:i
            beamBasis(i,:) = beamBasis(i,:).*(xvals-j*L/i*ones(1,n+2));
        end

    % Find current derivatives
    DbeamBasis(i,:) = diff(beamBasis(i,:))/deltax;
%   D2beamBasis(i,:) = diff(DbeamBasis(i,:))/deltax;

    % Recalculate stiffness matrix                      NOTE: CAN BE REWORKED TO BE MUCH MORE EFFICIENT
    for row = 1:i
        for col = 1:row
            product = DbeamBasis(row,:).*DbeamBasis(col,:);
            K(row,col) = trapz(product);
        end
    end
    K = K + K' - diag(diag(K));

    % Find condition number of current stiffness matrix
    [e_vecs, e_vals] = eig(K(1:i,1:i));
    e_vals = log(abs(diag(e_vals)));
    e_vals(e_vals == Inf | e_vals == -Inf) = 0;
    condnumbs_nodes(i) = max(e_vals) - min(e_vals);

end
disp('Condition numbers for polys found.')

% Find conditioning for Chebychev polys

K = zeros(shapes);
for i = 1:shapes

    % Chebyshev polysof the first kind
    beamBasis(i,:) = chebyshevT(2*i-1, xvals);

    % Find current derivatives
    DbeamBasis(i,:) = diff(beamBasis(i,:))/deltax;
%   D2beamBasis(i,:) = diff(DbeamBasis(i,:))/deltax;

    % Recalculate stiffness matrix                      NOTE: CAN BE REWORKED TO BE MUCH MORE EFFICIENT
    for row = 1:i
        for col = 1:row
            product = DbeamBasis(row,:).*DbeamBasis(col,:);
            K(row,col) = trapz(product);
        end
    end
    K = K + K' - diag(diag(K));

    % Find condition number of current stiffness matrix
    [e_vecs, e_vals] = eig(K(1:i,1:i));
    e_vals = log(abs(diag(e_vals)));
    e_vals(e_vals == Inf | e_vals == -Inf) = 0;
    condnumbs_cheby(i) = max(e_vals) - min(e_vals);

end
disp('Condition numbers for Chebys found.')

figure()
hold on;
plot(1:shapes, condnumbs_nodes, "LineStyle","none","Marker", ".")
plot(1:shapes, condnumbs_cheby, "LineStyle","none","Marker", ".")
title("Stiffness Matrix Condition Numbers vs. Shape Functions")
xlabel("# of Shape Functions")
ylabel("Condition Number")
legend("Polys with equally-spaced nodes", "Chebyshev polys of the first kind");
