function [rir_denoised_nm, edcs] = directional_denoise_SRIR(rir_noisy_nm,fs,pars,net)
%DIRECTIONAL_DENOISE_SRIR Reverberation tail de-noising for SRIRs.
%   rir_noisy_nm    : [numSmps, numSH]  - input in SHD, (N3D-ACN)
%   fs              : int               - Sampling Frequency
%   pars            : parameter struct
%   net             : DecayFitNet instance
[~,idx_max] = max(rms(rir_noisy_nm));
if idx_max == 1
    warning("SRIR might not be N3D normalized!")
end

N_sph = sqrt(length(rir_noisy_nm(1,:)))-1; % SH order
numSmps = size(rir_noisy_nm, 1);

if ~isfield(pars, 'fade')
    pars.fade = 'directional';  % or 'uniform' or sample ([1,numDirs])
end
if ~isfield(pars, 'numModes')
    pars.numModes = round((2^14)/size(rir_noisy_nm, 2));  % better numDirs
end
numBands = pars.numBands;

timeAxis = linspace(0, (numSmps - 1) / fs, numSmps );

numFadeSamples = round(0.25*fs);  % Cross-fade duration


%% Directional Decomposition of SRIR
[~, secDirs] = getTdesign(2*N_sph);
secPos = unitSph2cart(secDirs);
R = calculateRotationMatrix(secPos(1, :), [1, 0, 0]);
secPos = secPos * R;  % rotate to sec0 = [0,0]
[secAzi, secEle, secR] = cart2sph(secPos(:,1), secPos(:,2), secPos(:,3));
secDirs = [rad2deg(secAzi), rad2deg(secEle)];

numSecs = size(secDirs, 1);

%spatFilterCoeffs = sphButterworth(N_sph, 5, 2.5).';  % Or 'maxRE'
%spatFilterCoeffs = 'maxRE'
[A, B_AP] = designSphFilterBank(N_sph,secDirs,pars.spatFilterCoeffs,'AP');
[~, B_EP] = designSphFilterBank(N_sph,secDirs,pars.spatFilterCoeffs,'EP');


rirSecs = rir_noisy_nm * A;  % Beamformer Signals


%% DENOISE
rirSecsSynth = zeros(size(rirSecs));
n_max = zeros(1, numSecs);

measEdcSec = zeros(numSmps, numBands, numSecs);
estEdcSec = zeros(numSmps, numBands, numSecs);
synthEdcSec = zeros(numSmps, numBands, numSecs);


tic
for secIdx = 1:numSecs
    disp(["Direction " + num2str(secIdx)])
    rirChannel = rirSecs(:, secIdx);
    % estimate decayfit parameters
    % Last bool=true will also include the lowpass band below 125 Hz and the
    % highpass band above 8 kHz
    disp('=== Estimating EDC parameters with DecayFitNet ===');
    [t_est, a_est, n_est, norm_edc] = net.estimate_parameters(rirChannel,...
        true, true, pars.includeResidualBands); % bools: do_preprocess, do_scale_adjustment, includeResidualBands
    % scale A and N with norm_edc
    a_est = a_est .* norm_edc;
    n_est = n_est .* norm_edc;
    n_max(secIdx) = max(n_est);  % Store max over f

    [measEdc_, measNormEdc_] = rir2decay(rirChannel,...
        fs, net.filter_frequencies, true, true, true, pars.includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    
    
    % postprocessing: Treat bands with very little energy (<-120dB)
    t_est(a_est < 10e-13) = 0.01;
    t_difflast = (t_est(:, end) - t_est(:, end-1))./t_est(:, end-1);  %
    if any(abs(t_difflast) > 1)  % more than 100% change in HF probably noise
        t_est(abs(t_difflast) > 1, end) = t_est(abs(t_difflast) > 1, end-1);
    end


    estEdc_ = net.generate_synthetic_edcs(t_est, a_est, n_est, timeAxis).';
    
    measEdcSec(:, :, secIdx) = measEdc_ .* measNormEdc_;
    estEdcSec(:, :, secIdx) = estEdc_;

    % create a modal re-synthesis
    disp('Modal resynthesis:');
    estEdcClean = net.generate_synthetic_edcs(t_est, a_est, 0*n_est, timeAxis).';
    decayingNoises = modalSynthesisFromEdc(estEdcClean,net.filter_frequencies,pars.numModes,fs);
    modalRe = sum(decayingNoises,2);  % TODO move sum and sqrt

    disp('=== Getting EDC of modal resynthesis. ===');
    [modalReEdc_, modalReNorm_] = rir2decay(modalRe,...
        fs, net.filter_frequencies, true, true, true, pars.includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    synthEdcSec(:, :, secIdx) = modalReEdc_ .* modalReNorm_;
    
    % Write back
    rirSecsSynth(:,secIdx) = modalRe;
    
end
toc


%% Directional reconstruction of the denoised RIR

if strcmp(pars.fade, "directional") || strcmp(pars.fade, "uniform")
    % find cuton sample (t)
    % n_max needs to be scaled due to spatial decomposition
    n_max_total = 1/(trace(A'*A)/numSecs) * ... % orthonormalityerror
                    2*numBands*n_max;  % (simplified) bandpass scaling, 2b/fs
    
    smoothSecPow = medfilt1(rirSecs.^2, 100, [], 1);
    smoothSecPow = smoothSecPow ./ geometricSeries(n_max_total, numSmps);
    k_sec = zeros(1,numSecs);
    for secIdx = 1:numSecs
        k_sec(secIdx) = max([0.1*fs, ...
            find(smoothSecPow(:, secIdx) > (n_max_total(secIdx)*db2pow(6+12)), ...  % 12 for 200ms fade
                  1, 'last')]);
    end
    cutonSmp = k_sec;

    if strcmp(pars.fade, "uniform")
        cutonSmp = min(cutonSmp);  % same for all directions
    end
else
    cutonSmp = pars.fade;
end


%% create a random shift (negative delay), decorrelate, careful here
% shiftAmount = round(rand(1,numSecs)*0.000*fs);  % 0 deactivates
% rirSecsShift = zeros(size(rirSecs));
% for secIdx = 1:numSecs
%   rirSecsShift(:,secIdx) = rirSecs(:,secIdx);
%   rirSecsShift(1:shiftAmount(secIdx),secIdx) = 0;
%   rirSecsShift(:,secIdx) = circshift(rirSecsShift(:,secIdx), -shiftAmount(secIdx),1);
% end
% rirSecs = rirSecsShift;


%% Cross fade
ramp = linspace(0, 1, numFadeSamples);
fadeOut = sqrt(1-ramp);
fadeIn = sqrt(ramp);

rir_secs_faded = zeros(size(rirSecs));
rir_secs_faded(1:cutonSmp+numFadeSamples, :) = cat(1, rirSecs(1:cutonSmp, :), ...
     fadeOut.' .* rirSecs(cutonSmp+1: cutonSmp+numFadeSamples, :));

rir_secs_synth_faded = zeros(size(rirSecs));
rir_secs_synth_faded(cutonSmp+1: end, :) = cat(1, ...
    fadeIn.' .* rirSecsSynth(cutonSmp+1: cutonSmp+numFadeSamples, :), ...
    rirSecsSynth(cutonSmp+numFadeSamples+1 : end, :));


%% Should not be necessary: RMS norm
% if false
%     [~, k_dir] = max(abs(rir_in_nm(:, 1)));
%     % from peak+100ms to estimated noise floor sample
%     rms_norm = (rms(rir_secs_synth(k_dir+0.1*fs:cuton_smp, :)) ./ ...
%                       rms(rir_secs(k_dir+0.1*fs:cuton_smp, :)));
%     rir_secs_synth_faded = rir_secs_synth_faded ./ rms_norm;
% end


%% Transform back to SHD
% spat band limit
Ys = getSH(N_sph, deg2rad([secDirs(:,1), 90-secDirs(:,2)]), 'real').';
if strcmp(pars.spatFilterCoeffs,'maxRE')
    b_n = beamWeightsMaxEV(N_sph);
    c_n_an = b_n ./ sqrt((2*(0:N_sph).'+1) / (4*pi)); % remove m0 scaling
elseif strcmp(pars.spatFilterCoeffs, 'pwd')
    b_n = beamWeightsHypercardioid2Spherical(N_sph);
    c_n_an = b_n ./ sqrt((2*(0:N_sph).'+1) / (4*pi));  % remove m
else
    c_n_an = pars.spatFilterCoeffs;
end
Cspat = Ys.' * diag(replicatePerOrder(c_n_an)) * Ys;
assert(all(diag(Cspat) - ones(numSecs, 1)< 0.001, 'all'))
Cspat = Cspat ./ sum(Cspat, 1);

rir_secs_synth_faded = rir_secs_synth_faded * Cspat;

rir_denoised_nm = rir_secs_faded * B_AP' + ...
                  rir_secs_synth_faded * B_EP';

edcs = cat(4, measEdcSec, estEdcSec, synthEdcSec);


end