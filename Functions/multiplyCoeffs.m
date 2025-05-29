%% multiplyCoeffs.m
%  Joseph Anthony
%
% Created:         5/29/25
% Last Modified:   5/29/25
%
% Description: Multiplies two row polynomial coefficient vectors and
%   creates a new row coefficient vector
%
% INPUTS:
%   v1: coefficient vector 1
%   v2: coefficient vector 2
%   lengths: (optional) lengths of v1 and v2
% OUTPUTS:
%   prod: multiplied coefficient vector with lengths of v1 + v2

function product = multiplyCoeffs(v1, v2, vectorLengths)
    arguments 
        v1 (1,:) double
        v2 (1,:) double
        vectorLengths(1,2) double = [length(v1), length(v2)]
    end
    
    % Determine if a pad is needed
    pad = length(v1) + length(v2) - vectorLengths(1) - vectorLengths(2);

    % Determine which vector is shorter, if applicable, to reduce
    % iterations of for loop
    if vectorLengths(1) < vectorLengths(2)
        factorShort = v1;
        factorLong  = v2;
        lengthShort = vectorLengths(1);
    else
        factorShort = v2;
        factorLong  = v1;
        lengthShort = vectorLengths(2);
    end

    % Generate product vector
    product = zeros(1, vectorLengths(1) + vectorLengths(2));
    for i = 1:lengthShort
        product = product + [zeros(1,i-1),factorShort(i) * factorLong, zeros(1,lengthShort-i+1)];
    end

    % Pad with zeros as needed
    if pad ~= 0
        product = [product, zeros(1, pad)];
    end
end