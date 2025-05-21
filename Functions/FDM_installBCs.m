%% FDM_installBCs.m
%  Joseph Anthony
%
% Created:         5/21/25
% Last Modified:   5/21/25
%
% Description: Installs boundary conditions for a center FDM matrix with an
%   accuracy of order 2 and up to the fourth derivative
%
% INPUTS:
%   inputFDM:   input FDM matrix
%   order:      order of the input matrix
%   index:      index affected by the BCs; i = 0 or i = n+1 for endpoints
%   BCvector:   vector of BCs [u, Du, DDu, DDDu]
% OUTPUTS:
%   outputFDM:  FDM matrix with installed BCs

function outputFDM = FDM_installBCs(inputFDM, order, index, BCs)
    n = size(inputFDM, 2);
    
    end