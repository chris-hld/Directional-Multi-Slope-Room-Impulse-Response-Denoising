function [w_n] = sphButterworth(N_sph, k, l_c)
%SPHBUTTERWORTH Butterworth Filter, SH domain. Return modal weights w_n.
%   k : Filter order.
%   l_c : Cuton SH order.
% Chris Hold, 2022

w_n = 1./sqrt(1+([0:N_sph] ./ l_c).^(2*k));

% Unit Amplitude in Omega0
a = sum((2*[0:N_sph] + 1)/(4 * pi) .* w_n);
w_n = w_n / a;
end

