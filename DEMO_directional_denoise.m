close all; clear variables; clc; lastwarn('');

%% External dependencies
addpath(genpath('./DecayFitNet/'));
addpath(genpath('./SphFilterBank/'));
addpath(genpath('./matlabScripts'));

addpath(genpath('./Spherical-Harmonic-Transform/'));
addpath(genpath('./Spherical-Array-Processing/'));
addpath(genpath('./Higher-Order-Ambisonics/'));

%% Load RIR (ACN-N3D)
[rir_in_nm,fs_in] = audioread('./RIRs/eigenSRIR_doorway_6x10s.wav');  % SNR 61.3754

% MAKE SURE SRIRs ARE N3D!
warning("Converting Input to N3D")
rir_in_nm = convert_N3D_SN3D(rir_in_nm, 'sn2n');

%% Check version
versionchecker = DecayFitNetToolbox();
assert(strcmp(versionchecker.version, '0.0.7'), 'You are not using the correct DecayFitNet Toolbox version. Update submodule according to git commit.');
clear versionchecker

%% Parameters
nSlopes = 0;  % 0 for estimate nSlopes

pars = struct;

%pars.fade = 0.6 * pars.fs  % if you know where, you can specify the cuton sample

trim = 0.3;  % s (remove fade out, often present, makes estimation fail)

%% Preparing SRIR
pars.fs = fs_in;
assert(pars.fs==fs_in, 'fs of RIR does not match specified fs in parameters.');
N_sph_in = sqrt(length(rir_in_nm(1,:)))-1; % SH order

N_sph = 3;  % limit SH order if wanted
assert(N_sph<=N_sph_in, 'SH order exceeding input.');
rir_in_nm = rir_in_nm(:, 1:(N_sph+1)^2);  % Limit N_sph
rir_in_nm = rir_in_nm(1:end-trim*pars.fs,:); % remove fade out of RIR because it breaks the decay fit.
numSmps = size(rir_in_nm, 1);

rir_noisy_nm = rir_in_nm;

%% Setup Processing
pars.includeResidualBands = true;

pars.fBands = [125, 250, 500, 1000, 2000, 4000, 8000, 16000]; % set frequency bands for curve fitting
pars.numBands = numel(pars.fBands);
if pars.includeResidualBands; pars.numBands = pars.numBands+2; end
% get T,A values from network at octave bands
net = DecayFitNetToolbox(nSlopes, pars.fs);
net.filter_frequencies = pars.fBands;

pars.spatFilterCoeffs = sphButterworth(N_sph, 5, N_sph/2+1).';
%pars.spatFilterCoeffs = 'maxRE'

pars.secDirs = getSectorSteering(rir_noisy_nm, 'front');

%% Call
[rir_denoised_nm, edcs] = directional_denoise_SRIR(rir_noisy_nm,pars.fs,pars,net);

%% Metrics, and Plots
[SNR_diff,~,~] = compareMetrics(rir_noisy_nm,rir_denoised_nm,edcs,pars, 1);
[~,specDiffs,rt60diff] = compareMetrics(rir_in_nm,rir_denoised_nm,edcs,pars, 0);
disp("SNR improvement (dB) = " + SNR_diff)
disp("Mean Spectral Error (dB) = " + mean(mean(abs(specDiffs))));
disp("Mean rt60 Error (s) = " + mean(mean(abs(rt60diff))))

%% Save results
audiowrite("./RIRs/output/SRIRdenoisedN3D.wav", rir_denoised_nm, pars.fs,'BitsPerSample',32)
