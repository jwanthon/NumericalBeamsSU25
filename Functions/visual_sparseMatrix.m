%% visual_sparseMatrix.m
%  Joseph Anthony
%
% Created:         5/16/25
% Last Modified:   5/16/25
%
% Description: Creates a pre-formatted figure that plots the visuals for a
%   sparse square matrix.
%
% INPUTS:
%   matrix: numerical matrix (FEM, FDM, etc.)
% OUTPUTS:
%   fig:    figure


function fig = visual_sparseMatrix(matrix)
    % Create figure
    fig = figure;
    fig = heatmap(matrix);

    addpath('Data Files\');

    % Format and rescale colormap
    colormap = importdata('colormap_sparse.mat');
    fig.Colormap = colormap;
    [~, max] = bounds(matrix, "all");
    fig.ColorLimits = [-abs(max), abs(max)];

    fig.GridVisible = "off";
end