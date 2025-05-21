%% FDM_buildMtx.m
%  Joseph Anthony
%
% Created:         5/21/25
% Last Modified:   5/21/25
%
% Description: Builds a center FDM matrix with 2nd order accuracy.
%
% INPUTS:
%   d:      derivative
%   n:      mesh size
% OUTPUTS:
%   FDM:    FDM matrix

function FDM = FDM_buildMtx(derivative,meshsize)
    FDM = zeros(meshsize);
    switch derivative
        case 1
           v0 = zeros(1,meshsize);
           v1 = repmat(1/2,1,meshsize-1);
           FDM = FDM + diag(v0) + diag(v1,1) - diag(v1,-1);
        case 2
           v0 = repmat(-2,1,meshsize);
           v1 = ones(1,meshsize-1);
           FDM = FDM + diag(v0) + diag(v1,1) + diag(v1,-1);
        case 3
           v0 = zeros(1,meshsize);
           v1 = repmat(-1,1,meshsize-1);
           v2 = repmat(1/2,1,meshsize-2);
           FDM = FDM + diag(v0) + diag(v1,1) + diag(v2,2) ...
               - diag(v1, -1) - diag(v2, -2);
        case 4
           v0 = repmat(6,1,meshsize);
           v1 = repmat(-4,1,meshsize-1);
           v2 = ones(1,meshsize-2);
           FDM = FDM + diag(v0) + diag(v1,1) + diag(v2,2) ...
               + diag(v1, -1) + diag(v2, -2);
        case 5
           v0 = zeros(1,meshsize);
           v1 = repmat(5/2,1,meshsize-1);
           v2 = repmat(-2,1,meshsize-2);
           v3 = repmat(1/2,1,meshsize-3);
           FDM = FDM + diag(v0) + diag(v1,1) + diag(v2,2) + diag(v3, 3) ...
               - diag(v1, -1) - diag(v2, -2) - diag(v3, -3);
        case 6
           v0 = repmat(-20,1,meshsize);
           v1 = repmat(15,1,meshsize-1);
           v2 = repmat(-6,1,meshsize-2);
           v3 = ones(1,meshsize-3);
           FDM = FDM + diag(v0) + diag(v1,1) + diag(v2,2) + diag(v3, 3) ...
               + diag(v1, -1) + diag(v2, -2) + diag(v3, -3);
        otherwise
            warning('Input must be a scalar 1-6, not %d', derivative);
    end
end
