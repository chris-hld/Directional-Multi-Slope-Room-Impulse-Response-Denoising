function [G_filt] = tapering_compensation(w_n, f, N_target, soft_threshold, r)
%Frequency Gain of tapering compensation (diffuse field).
%   w_n : Modal weighting (tapering), column vector
%   f : Frequency Bins (Hz)
%   N_target :  Target SH order, e.g. 38
%   soft_threshold : soft clipping threshold in dB (max +6dB over)
%   r : Scatterer Radius, e.g. 0.0875
%
% Chris Hold, 2022
assert(iscolumn(w_n))

if nargin < 3
    N_target = 38;
end
if nargin < 4
    soft_threshold = 12;  % dB
end
if nargin < 5
    r = 0.0875;
end
c = 343;  % m/s
kr = 2*pi*f/c * r;

num_bins = length(kr);
N_tapered = length(w_n) - 1;
arrayType = 'rigid';
b_N_target = sphModalCoeffs(N_target, kr, arrayType);
b_N_tapered = sphModalCoeffs(N_tapered, kr, arrayType);


p_target = 1/(4*pi) * sqrt(sum(...
    repmat(2*(0:N_target)+1, num_bins, 1) .* abs(b_N_target).^2, 2));
p_tapered = 1/(4*pi) * sqrt(sum(repmat(w_n.', num_bins, 1) .* ...
    repmat(2*(0:N_tapered)+1, num_bins, 1) .* abs(b_N_tapered).^2, 2));


G_N = p_target ./ p_tapered;

% soft clip
G_filt = soft_clip(G_N, 10^(soft_threshold / 20)).';
end



function g_out = soft_clip(g_in, threshold)
%     Limit gain factor by soft clipping function. Limits gain factor to +6dB
%     beyond threshold point. (Pass values as factors/ratios, not dB!),
%     i.e. threshold = 10 ^ (dB / 20) .
%
%     Parameters
%     ----------
%     gain : array_like
%     threshold : float , e.g.(10 ^ (dB / 20))
%
%     Returns
%     -------
%     gain_clipped : array_like
gain = g_in / threshold;  % norm by threshold
gain(gain > 1) = 1 + tanh(gain(gain > 1) - 1);  % soft clipping to 2
g_out = gain * threshold;
end

