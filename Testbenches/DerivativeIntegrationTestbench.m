clc; clear all;
addpath('Functions\');
%% DerivativeIntegrationTestbench.m
%  Joseph Anthony
%
% Created:         6/20/25
% Last Modified:   6/20/25
%
% Description: Testbench for different integration and derivation techniques 
%   and their relative errors to determine what causes incorrect behaviors in
%   the FEM stiffness matrices
%
%% Diff vs. Gradient in K Matrix
clc; clear all;
n = 1000;         % Mesh size, number of interior points
shapes = 50;     % Number of shape functions, typically equal to n
L = 1;          % Beam length
modeCount = 5;  % Number of displayed modes

xvals       = linspace(0,L,n+2);
beamBasis   = zeros(shapes, n+2);
deltax      = L/(n+1);

% Find shape functions
for i = 1:shapes
    beamBasis(i,:) = xvals.*xvals.*(xvals - L * ones(1, n+2));
    for j = 1:i
        beamBasis(i,:) = beamBasis(i,:).*(xvals-j*L/i*ones(1,n+2));
    end
    beamBasis(i,:) = beamBasis(i,:)/max(abs(beamBasis(i,:)));   % Normalize shape function
end

% Find K using Gradient
DbeamBasis  = zeros(shapes,n+2);
D2beamBasis = zeros(shapes,n+2);

% Generate second derivatives
for i = 1:shapes
    DbeamBasis(i,:) = gradient(beamBasis(i,:),xvals);
    D2beamBasis(i,:) = gradient(DbeamBasis(i,:),xvals);
end

K1 = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = D2beamBasis(row,:).*D2beamBasis(col,:);
        K1(row,col) = trapz(product);
    end
end
K1 = K1 + K1' - diag(diag(K1));

% Find K using Diff
DbeamBasis  = zeros(shapes,n+1);
D2beamBasis = zeros(shapes,n);

% Generate second derivatives
for i = 1:shapes
    DbeamBasis(i,:) = diff(beamBasis(i,:))/deltax;
    D2beamBasis(i,:) = diff(DbeamBasis(i,:))/deltax;
end

K2 = zeros(shapes);
for row = 1:shapes
    for col = 1:row
        product = D2beamBasis(row,:).*D2beamBasis(col,:);
        K2(row,col) = trapz(product);
    end
end
K2 = K2 + K2' - diag(diag(K2));

% Error
error = abs(K1-K2);

fig2 = visual_sparseMatrix(error);
fig2.Title = "Absolute Error Between K1 and K2";

error1 = abs(error./K1)
fig3 = heatmap(error1);
fig3.Title = "Error Relative to K1";
fig3.GridVisible = 'off'

error2 = abs(error./K2)
fig4 = heatmap(error2);
fig4.Title = "Error Relative to K2";
fig4.GridVisible = 'off'


%% Diff vs. Gradient for Random Polys
clear all
hold off
L = 1;
xvals = linspace(-L/2,L/2,100);
max = 20;
n = 100;
deltax = xvals(2)-xvals(1);
% Generate random coefficients
coeff = ones(1,n);
for i = 1:n
    coeff(i) = randi(max)-max/2;
end

% Generate random exponents
for i = 1:n
    degrees(i) = randi(max);
end

% Generate y values
basis = zeros(1,length(xvals));
for i = 1:length(xvals)
    for j = 1:n
        basis(i) = basis(i) + coeff(j)*(xvals(i)^degrees(j));
    end
end
figure()
fig1 = plot(xvals,basis)
title(sprintf("Random Poly of Degree %d",n))

% Find derivative from gradient
grad_Dbasis = gradient(basis, xvals);
grad_D2basis = gradient(grad_Dbasis, xvals);
figure()
hold on
plot(xvals, grad_D2basis);

% Find derivative from diff
diff_Dbasis = diff(basis)/deltax;
diff_D2basis = diff(diff_Dbasis)/deltax;
plot(xvals(2:(length(xvals)-1)),diff_D2basis);

% Find exact derivative
exact_D2basis = 0*xvals;
for i = 1:length(xvals)
    for j = 1:n
        if degrees(j) == 2
            exact_D2basis(i) = exact_D2basis(i) + coeff(j)*2;
        elseif (degrees(j) == 1 || degrees(j) == 0)
            exact_D2basis(i) = exact_D2basis(i);
        else
            exact_D2basis(i) = exact_D2basis(i) + coeff(j)*degrees(j)*(degrees(j)-1)*(xvals(i)^(degrees(j)-2));
        end
    end
end
plot(xvals,exact_D2basis);
title(sprintf("Second Derivatives",n));
legend("gradient()", "diff()", "Exact");

% Find errors relative to exact
figure()
hold on
grad_error = grad_D2basis(2:end-1)./exact_D2basis(2:end-1);
diff_error = diff_D2basis./exact_D2basis(2:end-1);
plot(xvals(2:end-1), grad_error)
plot(xvals(2:end-1), diff_error)
title(sprintf("Relative Second Derivative Errors, 1 Trial",n));
legend("gradient()", "diff()")


%% Find accumulated errors for some set number of trials
trials = 1000;

grad_error = zeros(1, length(diff_D2basis));
diff_error = grad_error

for i = 1:trials
    % Generate random exponents
    for i = 1:n
        degrees(i) = randi(max);
    end
    % Generate y values
    basis = zeros(1,length(xvals));
    for i = 1:length(xvals)
        for j = 1:n
            basis(i) = basis(i) + coeff(j)*(xvals(i)^degrees(j));
        end
    end
    % Find derivative from gradient
    grad_Dbasis = gradient(basis, xvals);
    grad_D2basis = gradient(grad_Dbasis, xvals);
    % Find derivative from diff
    diff_Dbasis = diff(basis)/deltax;
    diff_D2basis = diff(diff_Dbasis)/deltax;
    exact_D2basis = 0*xvals;
    for i = 1:length(xvals)
        for j = 1:n
            if degrees(j) == 2
                exact_D2basis(i) = exact_D2basis(i) + coeff(j)*2;
            elseif (degrees(j) == 1 || degrees(j) == 0)
                exact_D2basis(i) = exact_D2basis(i);
            else
                exact_D2basis(i) = exact_D2basis(i) + coeff(j)*degrees(j)*(degrees(j)-1)*(xvals(i)^(degrees(j)-2));
            end
        end
    end
    grad_error = grad_error + grad_D2basis(2:end-1)./exact_D2basis(2:end-1);
    diff_error = diff_error + diff_D2basis./exact_D2basis(2:end-1);
end

% Plot average error
grad_error = grad_error/trials;
diff_error = diff_error/trials;
figure()
hold on;
plot(xvals(2:end-1), grad_error)
plot(xvals(2:end-1), diff_error)
title(sprintf("Average Relative Second Derivative Errors, %d Trials",trials));
legend("gradient()", "diff()");
%% Trapz vs. Direct Trap Sum Implementation for Random Polys
clear all
hold off
L = 1;
xvals = linspace(-L/2,L/2,100);
max = 20;
n = 100;
deltax = xvals(2)-xvals(1);
trials = 1000;

exact_vals = zeros(1, trials);
trapz_vals = exact_vals;
tsum_vals = exact_vals;

for iteration = 1:trials
    % Generate random coefficients
    coeff = ones(1,n);
    for i = 1:n
        coeff(i) = randi(max)-max/2;
    end
    
    % Generate random exponents
    for i = 1:n
        degrees(i) = randi(max);
    end
    
    % Generate y values
    basis = zeros(1,length(xvals));
    for i = 1:length(xvals)
        for j = 1:n
            basis(i) = basis(i) + coeff(j)*(xvals(i)^degrees(j));
        end
    end
    
    % Find trapz value
    trapz_int = trapz(xvals, basis);

    % Find implemented trapezoidal Riemann sum value
    tsum_int = deltax*(sum(basis)-1/2*basis(1)-1/2*basis(end));

    % Find exact integral value
    exact_int = 0;
    for i = 1:n
        exact_int = exact_int + coeff(i) * (xvals(end)^(degrees(i)+1)/(degrees(i)+1)-xvals(1)^(degrees(i)+1)/(degrees(i)+1));
    end

    exact_vals(iteration) = exact_int;
    trapz_vals(iteration) = trapz_int;
    tsum_vals(iteration) = tsum_int;

end

figure()
hold on
plot(1:trials, trapz_vals, "LineStyle","none", "Marker", ".", "Color", "red");
plot(1:trials, tsum_vals, "LineStyle","none", "Marker", "o", "Color", "blue");
plot(1:trials, exact_vals, "LineStyle","none", "Marker", "x", "Color", "black");
title("Integration Values for Random Polys");
xlabel("Trial #");
legend("trapz()", "Direct Trap Sum", "Exact Value");

figure()
hold on
trapz_err = trapz_vals./exact_vals;
tsum_err = tsum_vals./exact_vals;
plot(1:trials, trapz_err, "LineStyle","none", "Marker", ".", "Color", "red");
plot(1:trials, tsum_err, "LineStyle","none", "Marker", "o", "Color", "blue");
title("Relative Error for Different Integration Techniques");
xlabel("Trial #");
legend("trapz()", "Direct Trap Sum");
