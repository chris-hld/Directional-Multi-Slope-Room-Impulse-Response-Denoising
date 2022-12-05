function [s_out] = SHtoBin(sig_in_nm,fs,hrtf_path,TAPEREQ)
%SHtoBin - Renders real SH stream as binaural signals.
%
% Syntax:  [s_out_l,s_out_r] = SHtoBin(sig_in_nm,fs,hrtf_path)
%
% Inputs:
%    sig_in_nm  : [t x nSH] - real SH stream, ACN-N3D
%    fs         : int       - Sampling Frequency
%    hrtf_path  : string    - Path to HRTF set. Use THK_KU100 
%                             dataset for now.
%    TAPEREQ    : bool, opt - Taper and Diff-field EQ, default true
%
% Outputs:
%    s_out      : [t x 2]      - Binaural signals (l, r)
%
%
% Other m-files required: directSHT - Spherical-Harmonic-Transform 
%
% See also:
% ToDo: Check orthonormality limits, Check inverse factor
%
%   Chris Hold, 2022
%   Christoph.Hold@aalto.fi


if nargin<3
    hrtf_path = '~/data/HRTFs/THK_KU100/';
end
if nargin<4
    TAPEREQ = true;
end


addpath(hrtf_path)
load('HRIR_L2354.mat', 'HRIR_L2354')
%load HRIR_L2702.mat
load('KU100_CDFC.mat', 'hpcf')

nSH = size(sig_in_nm, 2);
order_SH = sqrt(nSH)-1;


set = HRIR_L2354;
assert(set.fs==fs)
azi = set.azimuth.';
colat = set.elevation.';
% quadweights sum to 1 instead of 4pi for this dataset
weights = 4*pi*set.quadWeight.';
ir_left = double(set.irChOne.');
ir_right = double(set.irChTwo.');
hrirTaps = size(ir_left, 2);

% transform hrir to SH domain
Y_N = getSH(order_SH, [azi, colat], 'real');
% perform SH transform
IR_nm_left = Y_N' * diag(weights) * ir_left;
IR_nm_right = Y_N' * diag(weights) * ir_right;

if TAPEREQ
    % Tapering
    numTapered = floor((order_SH-1)/2);
    c_taper = 0.5*(1-cos(2*pi * ...
        (numTapered+2:(2*numTapered+1))./(2*numTapered+2)));
    w_taper = [ones(order_SH-numTapered+1,1); c_taper.'];
    assert(length(w_taper) == order_SH+1)
    % Truncation EQ
    numEQTaps = 128;
    g_eq = tapering_compensation(w_taper,...
                                 linspace(0, fs/2, numEQTaps));
    b_truncEQ = fir2(numEQTaps, linspace(0, 1, numEQTaps), g_eq);
    [~, b_truncEQMin] = rceps(b_truncEQ);
    b_truncEQMin = b_truncEQMin(1:numEQTaps);  % Check tapering if necessary
    
    % apply SH Tapering
    sig_in_nm = replicatePerOrder(w_taper).' .* sig_in_nm;
    % zero pad and apply truncation EQ
    sig_in_nm = fftfilt(b_truncEQMin, [sig_in_nm; zeros(numEQTaps-1, nSH)]);
end

% zero pad input
sig_in_nm = [sig_in_nm; zeros(hrirTaps-1, nSH)];
% apply hrirs
s_out_l = sum(fftfilt(IR_nm_left.', sig_in_nm), 2);
s_out_r = sum(fftfilt(IR_nm_right.', sig_in_nm), 2);

% Headphone compensation
hpc_ir = hpcf.minPhase;
s_out_l = [s_out_l; zeros(length(hpc_ir) - 1, 1)];
s_out_r = [s_out_r; zeros(length(hpc_ir) - 1, 1)];

s_out_l = fftfilt(hpc_ir, s_out_l);
s_out_r = fftfilt(hpc_ir, s_out_r);

s_out = [s_out_l, s_out_r];
end

