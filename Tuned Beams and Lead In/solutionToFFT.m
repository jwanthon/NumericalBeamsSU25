%% solutionToAudio.m
%  Joseph Anthony
%
% Created:         5/8/25
% Last Modified:   5/12/25
%
% Description: Takes a solution vibration matrix and outputs an audio file 
%   representing the solution.
%
% INPUTS:
%   solution: solution vibration matrix
%   row:      chooses a row to compute the audio file from
%   timestep: time step across the horizontal axis
%   offset: amount of data points to discard on the leftmost side
% OUTPUTS:
%   fft:  FFT vector

function [fft_out, freqspace] = solutionToFFT(solution, row, timestep)

    samples = size(solution, 2);   % # of samples
    Fs = 1/timestep;                        % sampling frequency

    fft_out = fft(solution(row,:));
    freqspace = Fs/samples*(0:(samples-1));

end