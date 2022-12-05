function [A,B] = designSpatFilterBank(N,secDirs,pattern,mode)
%DESIGNSPHFILTERBANK SH domain Spatial Filter Bank, with amplitude or
%energy preservation during reconstruction.
% secDirs: [J, 2] - azi ele in deg
% pattern: 'maxRE', 'cardioid', 'pwd', or c_n coefficients
% mode: 'AP', 'EP'
%   apply with
%   s_sec = A * s_nm;
%   s_nm_out = B' * s_sec;
%
% Chris Hold, 2022

numSecs = size(secDirs, 1);
%Ys = getRSH(N, secDirs);
Ys = getSH(N, deg2rad([secDirs(:,1), 90-secDirs(:,2)]), 'real').';

if isstring(pattern) || ischar(pattern)
    switch pattern
        case 'cardioid'
            b_n = beamWeightsCardioid2Spherical(N);
            c_n = b_n ./ sqrt((2*(0:N).'+1) / (4*pi));  % remove m
        case 'maxRE'
            b_n = beamWeightsMaxEV(N);
            c_n = b_n ./ sqrt((2*(0:N).'+1) / (4*pi));  % remove m
        case 'pwd'
            b_n = beamWeightsHypercardioid2Spherical(N);
            c_n = b_n ./ sqrt((2*(0:N).'+1) / (4*pi));  % remove m
        otherwise
            error("Pattern not implemented");
    end
else
    c_n = pattern;
end
assert(all(size(c_n) == [N+1, 1]))

% Construct Analysis matrix
A = diag(replicatePerOrder(c_n)) * Ys;

% Construct Synthesis matrix
B = [];
switch mode
    case 'AP'
        beta = (sqrt(4*pi)) / (A(1, 1) * numSecs);
        B = beta * diag(replicatePerOrder(1./(c_n/c_n(1)))) * Ys;
    case 'EP'
        beta =  (4*pi) / (A(:,1)' * A(:,1) * numSecs);
        B = sqrt(beta) * diag(replicatePerOrder(1./(c_n/c_n(1)))) * Ys;
    otherwise
        warning("Not implemented")
end


end

