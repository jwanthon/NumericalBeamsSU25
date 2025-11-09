clc; clear all;
addpath('Functions\');

%% StringExamplesFEM.m
%  Joseph Anthony
%
% Created:         11/08/25
% Last Modified:   11/08/25
%
% Description: Produces several diagrams of DD and DN string mode shapes
%   and basis functions.

n = 100;       % Mesh size, number of interior points
shapes = n;   % Number of shape functions, typically equal to n
L = 1;          % Beam length
modecount = 5;  % Number of displayed modes

%% Both DD, Analytic D-N
xvals = linspace(0,L,n+2);
deltax = xvals(2)-xvals(1);
Dxvals = zeros(1, n+1); % Generate a mesh of xvals of derivative
for i = 1:length(xvals) - 1
    Dxvals(i) = (xvals(i) + xvals(i+1))/2;
end

% Initialize variables
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes,n+1);
basis_DDmodes = basis;
basis_DNmodes = basis;
basis_tents = basis;
basis_waves = basis;
Dbasis_DDmodes = Dbasis;
Dbasis_DNmodes = Dbasis;
Dbasis_tents = Dbasis;
Dbasis_waves = Dbasis;

% Generate shape functions and their respective derivatives
for i = 1:shapes
    basis_DDmodes(i,:) = sin(pi*i*xvals);
    Dbasis_DDmodes(i,:) = diff(basis_DDmodes(i,:))/deltax;
    
    basis_DNmodes(i,:) = sin((2*i-1)/2*pi*xvals);
    Dbasis_DNmodes(i,:) = diff(basis_DNmodes(i,:))/deltax;

    basis_tents(i,i+1) = 1;
    Dbasis_tents(i,:) = diff(basis_tents(i,:))/deltax;
end

% Generate stiffness matrices
K_DDmodes = zeros(shapes);
K_DNmodes = zeros(shapes);
K_tents = zeros(shapes);

for row = 1:shapes
    for col = 1:row
        K_DDmodes(row,col) = trapz(Dxvals, Dbasis_DDmodes(row,:).*Dbasis_DDmodes(col,:));
        K_DNmodes(row,col) = trapz(Dxvals, Dbasis_DNmodes(row,:).*Dbasis_DNmodes(col,:));
        K_tents(row,col) = trapz(Dxvals, Dbasis_tents(row,:).*Dbasis_tents(col,:));
    end
end

K_DDmodes = K_DDmodes + K_DDmodes' - diag(diag(K_DDmodes));
K_DNmodes = K_DNmodes + K_DNmodes' - diag(diag(K_DNmodes));
K_tents = K_tents + K_tents' - diag(diag(K_tents));

% Determine sorted eigenbasis
[evecs_DDmodes, evals_DDmodes] = eig(K_DDmodes);
[~, index] = sort(diag(abs(evals_DDmodes))); 
evals_DDmodes = evals_DDmodes(index,index);
evecs_DDmodes = evecs_DDmodes(:, index);

[evecs_DNmodes, evals_DNmodes] = eig(K_DNmodes);
[~, index] = sort(diag(abs(evals_DNmodes))); 
evals_DNmodes = evals_DNmodes(index,index);
evecs_DNmodes = evecs_DNmodes(:, index);

[evecs_tents, evals_tents] = eig(K_tents);
[~, index] = sort(diag(abs(evals_tents))); 
evals_tents = evals_tents(index,index);
evecs_tents = evecs_tents(:, index);



% Compute mode shapes
modes_DDmodes = basis_DDmodes'*evecs_DDmodes;
modes_DDmodes = modes_DDmodes';

modes_DNmodes = basis_DNmodes'*evecs_DNmodes;
modes_DNmodes = modes_DNmodes';

modes_tents = basis_tents'*evecs_tents;
modes_tents = modes_tents';

figure();
hold on
for i = 1:modecount
    modes_tents(i,:) = modes_tents(i,:)/max(abs(modes_tents(i,:)));
    plot(xvals, modes_tents(i,:));
end
xlabel("String Position");
ylabel("Relative Displacement");
legend('1','2','3','4','5');


figure();
hold on
for i = 1:modecount
    modes_DDmodes(i,:) = modes_DDmodes(i,:)/max(abs(modes_DDmodes(i,:)));
    plot(xvals, modes_DDmodes(i,:));
end
xlabel("String Position");
ylabel("Relative Displacement");
legend('1','2','3','4','5');


figure();
hold on
for i = 1:modecount
    modes_DNmodes(i,:) = modes_DNmodes(i,:)/max(abs(modes_DNmodes(i,:)));
    plot(xvals, modes_DNmodes(i,:));
end
xlabel("String Position");
ylabel("Relative Displacement");
legend('1','2','3','4','5');


%% Numerical D-N
dnshapes = shapes+1;

% Initialize variables
basis_numdn = zeros(dnshapes, n+2);
Dbasis_numdn = zeros(dnshapes,n+1);


% Generate shape functions and their respective derivatives
for i = 1:dnshapes
    basis_numdn(i,i+1) = 1;
    Dbasis_numdn(i,:) = diff(basis_numdn(i,:))/deltax;
end

% Generate stiffness matrices
K_numdn = zeros(dnshapes);

for row = 1:dnshapes
    for col = 1:row
        K_numdn(row,col) = trapz(Dxvals, Dbasis_numdn(row,:).*Dbasis_numdn(col,:));
    end
end

K_numdn = K_numdn + K_numdn' - diag(diag(K_numdn));

% Determine sorted eigenbasis
[evecs_numdn, evals_numdn] = eig(K_numdn);
[~, index] = sort(diag(abs(evals_numdn))); 
evals_numdn = evals_numdn(index,index);
evecs_numdn = evecs_numdn(:, index);

% Compute mode shapes
modes_numdn = basis_numdn'*evecs_numdn;
modes_numdn = modes_numdn';

figure();
hold on
modes_numdn(1,:)=-modes_numdn(1,:);
for i = 1:modecount
    modes_numdn(i,:) = modes_numdn(i,:)/max(abs(modes_numdn(i,:)));
    plot(xvals, modes_numdn(i,:));
end
xlabel("String Position");
ylabel("Relative Displacement");
legend('1','2','3','4','5');

%% Plot eigenvalue spectrum
espectrum_DDmodes = abs(diag(evals_DDmodes))/abs(evals_DDmodes(1,1));
rootvals_DDmodes = sqrt(espectrum_DDmodes);

espectrum_tents = abs(diag(evals_tents))/abs(evals_tents(1,1));
rootvals_tents = sqrt(espectrum_tents);

espectrum_DNmodes = abs(diag(evals_DNmodes))/abs(evals_DNmodes(1,1));
rootvals_DNmodes = sqrt(espectrum_DNmodes);

espectrum_numdnmodes = abs(diag(evals_numdn))/abs(evals_numdn(1,1));
rootvals_numdnmodes = sqrt(espectrum_numdnmodes);

figure();
hold on
plot(rootvals_DDmodes, "LineStyle", "none", "Marker", "o");
plot(rootvals_tents, "LineStyle", "none", "Marker", "*");
xlim([1 50]);
legend("Analytic D-D Modes", "Tent D-D Modes", 'Location', 'southeast')


figure();
hold on
plot(rootvals_DNmodes, "LineStyle", "none", "Marker", "o");
plot(rootvals_numdnmodes, "LineStyle", "none", "Marker", "*");
xlim([1 50]);
legend("Analytic D-N Modes", "Tent D-N Modes", 'Location', 'southeast')