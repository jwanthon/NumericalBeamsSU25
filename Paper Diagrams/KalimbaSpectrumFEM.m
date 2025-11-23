clc; clear all;
addpath('Functions\');

%% KalimbaSpectrumFEM.m
%  Joseph Anthony
%
% Created:         11/23/25
% Last Modified:   11/23/25
%
% Description: Creates frequency lines for the whammy bar kalimba
%% Parameters
n = 100;            % Mesh size, number of interior points
shapes = n;         % Number of shape functions

L = 1;              % Beam length
modeCount = 5;      % Number of displayed modes

lowfreqs = 5;       % Number of frequency lines
resolution = 100;   % Number of points

%% Multi-Mode

xvals = linspace(0,L,n+2);
deltax = mean(diff(xvals));
basis = zeros(shapes, n+2);
Dbasis = zeros(shapes, n+1);
D2basis = zeros(shapes, n);

freqlines = zeros(lowfreqs, resolution);

bvals = linspace(L/8, 3*L/8, resolution);

for ksample = 1:resolution
    
    b = bvals(ksample);

    % Generate shape functions and their respective derivatives
    for i = 1:shapes
        basis(i,:) = xvals.*(xvals-b*ones(1,n+2)).*sin((2*i-1)*pi/2*xvals);
        Dbasis(i,:) = diff(basis(i,:))/deltax;
        D2basis(i,:) = diff(Dbasis(i,:))/deltax;
    end
    
    % Generate symmetric stiffness matrices
    K = zeros(shapes);
    M = zeros(shapes);
    
    for row = 1:shapes
        for col = 1:row
            K(row,col) = trapz(xvals(2:end-1), D2basis(row,:).*D2basis(col,:));
        end
    end
    K = K + K' - diag(diag(K));
    
    % Determine sorted eigenbasis
    [evecs, evals] = eig(K);
    [~, index] = sort(diag(abs(evals))); 
    evals = diag(abs(evals(index,index)));

    freqlines(:,ksample) = evals(1:lowfreqs);

end

hold on
grid on
for i = 1:lowfreqs
    plot(bvals,log(freqlines(i,:)));
end
legend('1','2','3','4','5', 'Location', 'southeast');
