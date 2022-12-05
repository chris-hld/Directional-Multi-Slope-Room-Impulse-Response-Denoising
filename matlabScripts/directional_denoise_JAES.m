close all; clear variables; clc; lastwarn('');

%% Flags
PLOT = true;
BIN = true;

%% External dependencies
addpath(genpath('../DecayFitNet/'));
addpath(genpath('../SphFilterBank/'));

addpath(genpath('../RIRs'));

%% Set Paths
% Dependencies
gitpath = '~/git';

libs = {'Higher-Order-Ambisonics/',...
    'Spherical-Harmonic-Transform/',...
    'Spherical-Array-Processing/',...
    'Array-Response-Simulator/'};

for lib = libs
    addpath(fullfile(gitpath, lib{1}));
end

%% Check version
versionchecker = DecayFitNetToolbox();
assert(strcmp(versionchecker.version, '0.0.7'), 'You are not using the correct DecayFitNet Toolbox version. Update submodule according to git commit.');
clear versionchecker

%% Parameters
fs = 48000;

target_snr = 300; % dB
N_sph = 3;

numModes = 2^13;
includeResidualBands = true;
nSlopes = 0;  % 0 for estimate nSlopes

%% Load RIR
% Bug: eigenSRIR was not n3d
%[rir_in_nm,fs_in] = audioread('eigenSRIR_doorway_1x60s.wav');  % SNR 60.5091
[rir_in_nm,fs_in] = audioread('eigenSRIR_doorway_6x10s.wav');  % SNR 61.3754, this seems ok
%[rir_in_nm,fs_in] = audioread('eigenSRIR_doorway_60x1s.wav');  % SNR61.8439, don't use this
%[rir_in_nm,fs_in] = audioread('RIR_250cm.wav');  % SNR 46.3595
%[rir_in_nm,fs_in] = audioread('shoeboxHybridRIR-SH5.wav');

rir_in_nm = rir_in_nm ./ (max(rir_in_nm(:, 1)) *2);

trim = 0.3;  % s (remove fade out)

assert(fs==fs_in, 'fs of RIR does not match specified fs in parameters.');
N_sph_in = sqrt(length(rir_in_nm(1,:)))-1; % SH order
assert(N_sph<=N_sph_in, 'SH order exceeding input.');
rir_in_nm = rir_in_nm(:, 1:(N_sph+1)^2);  % Limit N_sph
rir_in_nm = rir_in_nm(1:end-trim*fs,:); % remove fade out of RIR because it breaks the decay fit.
numSmps = size(rir_in_nm, 1);

% Load noise
noise_nm = randn(size(rir_in_nm)); % create noise, diffuse band-limited
%noise_nm = audioread("EM_backgroundNoise.wav");
%noise_nm = noise_nm(fs:fs+numSmps-1, 1:(N_sph+1)^2);


%% Setup
time = (1:numSmps)';

fbands = [125, 250, 500, 1000, 2000, 4000, 8000, 16000]; % set frequency bands for curve fitting
numBands = numel(fbands);
if includeResidualBands; numBands = numBands+2; end

% get T,A values from network at octave bands
net = DecayFitNetToolbox(nSlopes, fs);

net.filter_frequencies = fbands;

time_axis = linspace(0, (numSmps - 1) / fs, numSmps );

%% add noise to RIR
noise_energy = 1/(((N_sph+1)^2) *...
                length(noise_nm))*trace((noise_nm.' * noise_nm));
% TODO: RMS
% to define signal-to-noise ratio
signal_energy = 1/(((N_sph+1)^2)*0.5*fs)*...
                 trace(rir_in_nm(1:0.5*fs,:).' * rir_in_nm(1:0.5*fs,:));
% define noise accordingly
noise_gain =  sqrt(signal_energy)/sqrt(noise_energy) * 10^(-target_snr/20); % a larger number means a noisier RIR, or choose 0 for no noise
noise_nm = (noise_nm * noise_gain);
noise_energy_scaled = 1/(((N_sph+1)^2) *...
                length(noise_nm))*trace((noise_nm.' * noise_nm));

rir_noisy_nm = rir_in_nm + noise_nm; % add in noise here

fprintf('Signal energy: %.02f dB, Add. noise energy scaled: %.02f dB. (SNR: %.02f dB).\n', ...
        pow2db(signal_energy), pow2db(noise_energy_scaled), target_snr);

%% Directional Decomposition of SRIR
[~, secDirs] = getTdesign(2*N_sph);
secPos = unitSph2cart(secDirs);
R = calculateRotationMatrix(secPos(1, :), [1, 0, 0]);
secPos = secPos * R;  % rotate to sec0 = [0,0]
[secAzi, secEle, secR] = cart2sph(secPos(:,1), secPos(:,2), secPos(:,3));
secDirs = [rad2deg(secAzi), rad2deg(secEle)];

numSecs = size(secDirs, 1);

spatFilterCoeffs = sphButterworth(N_sph, 5, 2.5).';  % Or 'maxRE'
%spatFilterCoeffs = 'maxRE'
[A, B_AP] = designSphFilterBank(N_sph,secDirs,spatFilterCoeffs,'AP');
[~, B_EP] = designSphFilterBank(N_sph,secDirs,spatFilterCoeffs,'EP');


rir_secs = rir_noisy_nm * A;  % Sector Signals

[max_secamp, max_secidx] =  max(max(abs(rir_secs))); % check if the amplitude is highest at the directions where the source is

%% DENOISE

rir_secs_synth = zeros(size(rir_secs));
switch_sample_dir = zeros(numSecs, 1);
n_max = [];

meas_edc_sec = zeros(numSmps, numBands, numSecs);
est_edc_sec = zeros(numSmps, numBands, numSecs);
synth_edc_sec = zeros(numSmps, numBands, numSecs);


tic
for secIdx = 1:numSecs
    disp(["Sector " + num2str(secIdx)])
    rir_channel = rir_secs(:, secIdx);
    % estimate decayfit parameters
    % Last bool=true will also include the lowpass band below 125 Hz and the
    % highpass band above 8 kHz
    disp('=== Estimating EDC parameters with DecayFitNet ===');
    [t_values, a_values, n_values, norm_edc] = net.estimate_parameters(rir_channel, true, true, includeResidualBands); % bools: do_preprocess, do_scale_adjustment, includeResidualBands

    [meas_edc_, meas_norm_edc_] = rir2decay(rir_channel, fs, net.filter_frequencies, true, true, true, includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    
    % scale A and N with norm_edc
    a_values = a_values .* norm_edc;
    n_values = n_values .* norm_edc;

    % postprocessing: Treat bands with very little energy (<-120dB)
    t_values(a_values < 10e-13) = 0.01;
    t_diff_last = (t_values(:, end) - t_values(:, end-1))./t_values(:, end-1);  %
    if any(abs(t_diff_last) > 1)  % more than 100% change in HF probably noise
        t_values(abs(t_diff_last) > 1, end) = t_values(abs(t_diff_last) > 1, end-1);
    end

    n_max = max([n_values, n_max]);  % TODO, store per sector

    est_edc_ = net.generate_synthetic_edcs(t_values, a_values, n_values, time_axis).';
    
    meas_edc_sec(:, :, secIdx) = meas_edc_ .* meas_norm_edc_;
    est_edc_sec(:, :, secIdx) = est_edc_;

    % create a modal re-synthesis
    disp('Modal resynthesis:');
    est_edc_no_noise = net.generate_synthetic_edcs(t_values, a_values, 0*n_values, time_axis).';
    decayingNoises = modalSynthesisFromEdc(est_edc_no_noise,net.filter_frequencies,numModes,fs);
    modalRe = sum(decayingNoises,2);  % TODO move sum and sqrt

    disp('=== Getting EDC of modal resynthesis. ===');
    [modalRe_edc_, modalRe_norm_edc] = rir2decay(modalRe, fs, net.filter_frequencies, true, true, true, includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    synth_edc_sec(:, :, secIdx) = modalRe_edc_ .* modalRe_norm_edc;
    
    % Write back
    rir_secs_synth(:,secIdx) = modalRe;
    
end
toc


%% Directional reconstruction of the denoised RIR
% find cuton sample (t)
% n_max needs to be scaled due decomposition
n_max_total = 1/(trace(A'*A)/numSecs) * ... % orthonormalityerror
                2*numBands*n_max;  % (simplified) bandpass scaling, 2b/fs
%n_max_total = 2*numBands*n_max
%n_max_total = 2 * n_max_total;  % TODO
smoothSecPow = medfilt1(rir_secs.^2, 100, [], 1);
smoothSecPow = smoothSecPow ./ geometricSeries(n_max_total, numSmps);
for secIdx = 1:numSecs
    k_sec(secIdx) = max([0.1*fs, find(smoothSecPow(:, secIdx) > (n_max_total*db2pow(6+12)), ...  % 12 for 200ms fade
              1, 'last')]);
end
cuton_smp = k_sec;
cuton_smp = max(k_sec);  % same for all sectors


%% create a random shift (negative delay), decorrelate, careful here
shiftAmount = round(rand(1,numSecs)*0.000*fs);
rir_secs_shift = zeros(size(rir_secs));
for secIdx = 1:numSecs
  rir_secs_shift(:,secIdx) = rir_secs(:,secIdx);
  rir_secs_shift(1:shiftAmount(secIdx),secIdx) = 0;
  rir_secs_shift(:,secIdx) = circshift(rir_secs_shift(:,secIdx), -shiftAmount(secIdx),1);
end


%% Cross fade
num_fade_samples = round(0.3*fs);  % Cross-fade duration
ramp = linspace(0, 1, num_fade_samples);
fade_out = sqrt(1-ramp);
fade_in = sqrt(ramp);


rir_secs_faded = zeros(size(rir_secs_shift));
rir_secs_faded(1:cuton_smp+num_fade_samples, :) = cat(1, rir_secs_shift(1:cuton_smp, :), ...
     fade_out.' .* rir_secs_shift(cuton_smp+1: cuton_smp+num_fade_samples, :));

rir_secs_synth_faded = zeros(size(rir_secs_shift));
rir_secs_synth_faded(cuton_smp+1: end, :) = cat(1, ...
    fade_in.' .* rir_secs_synth(cuton_smp+1: cuton_smp+num_fade_samples, :), ...
    rir_secs_synth(cuton_smp+num_fade_samples+1 : end, :));


% only for comparison
rir_nm_faded = zeros(size(rir_noisy_nm));
rir_nm_faded(1:cuton_smp+num_fade_samples, :) = cat(1, rir_noisy_nm(1:cuton_smp, :), ...
     fade_out.' .* rir_noisy_nm(cuton_smp+1: cuton_smp+num_fade_samples, :));


%% Should not be necessary: RMS sector norm
if false
    [~, k_dir] = max(abs(rir_in_nm(:, 1)));
    % from peak+100ms to estimated noise floor sample
    rms_norm = (rms(rir_secs_synth(k_dir+0.1*fs:cuton_smp, :)) ./ ...
                      rms(rir_secs(k_dir+0.1*fs:cuton_smp, :)));
    rir_secs_synth_faded = rir_secs_synth_faded ./ rms_norm;
end


%% Transform back to SHD
% spat band limit
Ys = getRSH(N_sph, secDirs);
if strcmp(spatFilterCoeffs,'maxRE')
    b_n = beamWeightsMaxEV(N_sph);
    c_n_an = b_n ./ sqrt((2*(0:N_sph).'+1) / (4*pi)); % remove m0 scaling
elseif strcmp(spatFilterCoeffs, 'pwd')
    b_n = beamWeightsHypercardioid2Spherical(N_sph);
    c_n_an = b_n ./ sqrt((2*(0:N_sph).'+1) / (4*pi));  % remove m
else
    c_n_an = spatFilterCoeffs;
end
Cspat = Ys.' * diag(replicatePerOrder(c_n_an)) * Ys;
assert(all(diag(Cspat) - ones(numSecs, 1)< 0.001, 'all'))
Cspat = Cspat ./ sum(Cspat, 1);

rir_secs_synth_faded = rir_secs_synth_faded * Cspat;

rir_denoised_nm = rir_secs_faded * B_AP' + ...
                  rir_secs_synth_faded * B_EP';

% for comparison (not part of processing)
rir_synth_nm_ = rir_secs_synth * Cspat * B_EP';
rir_secs_in_ = rir_in_nm * A;
rir_secs_out_ = rir_denoised_nm * A;
%rir_denoised_nm = rir_nm_faded + ...
%    rir_secs_synth_faded * B_EP';


%% Metrics
lateIdx = round(0.1 * fs);

[~, k_max] = max(abs(rir_in_nm(:, 1)));
k_max = k_max - 10;  % get the first couple of onset samples
k_range = round(0.1 * fs);

%[meas_edc0] = rir2decay(sqrt(4*pi)*rir_in_nm(:,1), fs, fbands, true, true, false, includeResidualBands);
%[synth_edc0] = rir2decay(sqrt(4*pi)*rir_synth_nm(:,1), fs, fbands, true, true, false, includeResidualBands);

secRT60 = zeros(1, numSecs);
measRT60sec = zeros(numSecs, numBands);
estRT60sec = zeros(numSecs, numBands);
for idxSec = 1:numSecs
    meas_edc_sec_ = meas_edc_sec(:, :, idxSec);
    meas_edc_sec_ = meas_edc_sec_ ./ max(meas_edc_sec_, [], 1);
    est_edc_sec_ = est_edc_sec(:, :, idxSec);
    est_edc_sec_ = est_edc_sec_ ./ max(est_edc_sec_, [], 1);
    for idxBand = 1:numBands
        measRT60sec(idxSec,idxBand) = 2*(find(meas_edc_sec_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
        estRT60sec(idxSec,idxBand) = 2*(find(est_edc_sec_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
    end
end


meas_edc_sec_in = zeros(size(meas_edc_sec));
meas_edc_sec_out = zeros(size(meas_edc_sec));
measRT60secIn = zeros(numSecs, numBands);
measRT60secOut = zeros(numSecs, numBands);
for idxSec = 1:numSecs
    [meas_edc_sec_in_, meas_norm_edc_sec_in_] = rir2decay(rir_secs_in_(:,idxSec), fs, net.filter_frequencies, true, true, true, includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    [meas_edc_sec_out_, meas_norm_edc_sec_out_] = rir2decay(rir_secs_out_(:,idxSec), fs, net.filter_frequencies, true, true, true, includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    meas_edc_sec_in(:, :, idxSec) = meas_edc_sec_in_ .* meas_norm_edc_sec_in_;
    meas_edc_sec_out(:, :, idxSec) = meas_edc_sec_out_ .* meas_norm_edc_sec_out_;
    
    for idxBand = 1:numBands
        if meas_norm_edc_sec_in_(idxBand) < 10e-13  % postpro: no energy in band
            measRT60secIn(idxSec,idxBand) = 0.01;
            measRT60secOut(idxSec,idxBand) = 0.01;
        else
            measRT60secIn(idxSec,idxBand) = 2*(find(meas_edc_sec_in_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
            measRT60secOut(idxSec,idxBand) = 2*(find(meas_edc_sec_out_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
        end
    end
end

for idxSec=1:numSecs
[h_secs_in(:,idxSec),f] = freqz(rir_secs_in_(lateIdx:end,idxSec),1,2^12,fs);
[h_secs_denoised(:,idxSec),~] = freqz(rir_secs_faded(lateIdx:end,idxSec) + ...
                           rir_secs_synth_faded(lateIdx:end,idxSec),...
                           1,2^12,fs);
end

% 4 pi might be missing here, but should cancel
SNR_in = trace(rir_noisy_nm(k_max:k_max+k_range, :)'*rir_noisy_nm(k_max:k_max+k_range, :)) / ...
    trace(rir_noisy_nm(end-k_range:end, :)'*rir_noisy_nm(end-k_range:end, :));
SNR_out = trace(rir_denoised_nm(k_max:k_max+k_range, :)'*rir_denoised_nm(k_max:k_max+k_range, :)) / ...
    trace(rir_denoised_nm(end-k_range:end, :)'*rir_denoised_nm(end-k_range:end, :));
SNR_improvement = 10*log10(SNR_out/SNR_in);

for idxSec=1:numSecs
    specErrors(idxSec) = mean(abs(fracBandSmooth(mag2db(abs(h_secs_denoised(:,idxSec))), 3) - ...
        fracBandSmooth(mag2db(abs(h_secs_in(:,idxSec))), 3)));
end
disp("SNR improvement (dB) = " + SNR_improvement)
disp("Mean Spectral Error (dB) = " + mean(specErrors))
disp("Mean rt60 Error (s) = " + mean(mean(abs(measRT60secOut - measRT60secIn))))


%% Analyze and Plot
if PLOT
defFigsize = [3.25, 2.5];
plotBandIdx = 5;
plotSecIdx = [1, 8, 11, 15, 18, 23];
if includeResidualBands
    fbandsString = ["Lo" num2cell(fbands) "Hi"];
else
    fbandsString = string(num2cell(fbands));
end


%%% IR
figure;
hold on;
plot(time_axis, pow2db(dot(rir_noisy_nm, rir_noisy_nm, 2)), ':');
plot(time_axis, pow2db(dot(rir_denoised_nm, rir_denoised_nm, 2)), '-');
xline(cuton_smp/fs, '--', "Cross-Over")
yline(pow2db(n_max_total), ':', "Noise-Est")
ylim([-100, 0])
xlabel("Time (s)")
ylabel("Magnitude (dB)")
legend("Input", "Output")
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/sigs.eps", 'Resolution',800)


% figure;
% sgtitle('RIR p(t)')
% 
% subplot(3, 1, 1)
% plot(time_axis, sqrt(4*pi) * rir_in_nm(:,1));
% hold on;
% plot(time_axis, sqrt(4*pi) * rir_noisy_nm(:,1), '--');
% legend("clean RIR", "noisy RIR")
% ylim([-1, 1])
% subtitle("IN")
% 
% subplot(3, 1, 2)
% plot(time_axis, sqrt(4*pi) * rir_synth_nm(:,1));
% ylim([-1, 1])
% subtitle("OUT Synth")
% 
% subplot(3, 1, 3)
% plot(time_axis, sqrt(4*pi) * rir_denoised_nm(:,1));
% hold on
% xline(cuton_smp/fs, '-', "Cross-Over")
% ylim([-1, 1])
% subtitle("Synthesized")

figure;
sgtitle('RIR p(t)')

subplot(3, 1, 1)
plot(time_axis, db(sqrt(4*pi) * rir_in_nm(:,1)));
hold on;
plot(time_axis, db(sqrt(4*pi) * rir_noisy_nm(:,1)), ':');
ylim([-100, 0])
yline(pow2db(n_max_total), '-', "noise-est")
grid on
%legend("clean RIR", "noisy RIR")
subtitle("IN")

subplot(3, 1, 2)
plot(time_axis, db(sqrt(4*pi) * rir_synth_nm_(:,1)));
ylim([-100, 0])
grid on
subtitle("OUT Synth")

subplot(3, 1, 3)
plot(time_axis, db(sqrt(4*pi) * rir_denoised_nm(:,1)));
hold on
xline(cuton_smp/fs, '-', "Cross-Over")
xlabel("time in s")
ylabel("in dB")
ylim([-100, 0])
grid on
subtitle("Synthesized")


%%% EDC
% figure; 
% sgtitle("EDC p(t), meas (solid) vs. synth (dashed)")
% 
% plot(time_axis, pow2db(meas_edc0))
% hold on
% set(gca,'ColorOrderIndex',1)
% plot(time_axis, pow2db(synth_edc0), '--')
% legend(split(fbandsString))
% xlim([0, time_axis(end)])
% ylim([max(pow2db(meas_edc0),[],'all')-100, max(pow2db(meas_edc0),[],'all')])
% xlabel("time (s)")
% ylabel("in dB")
% grid on

figure;
hold on
%sec_norm_ = max(meas_edc_sec(:, :, 1), [], 1);  % norm on 1st sec
for idxSec = 1:length(plotSecIdx)
    meas_edc_sec_in_ = meas_edc_sec(:, :, plotSecIdx(idxSec));
    sec_norm_ = max(meas_edc_sec_in_, [], 1);
    meas_edc_sec_in_ = meas_edc_sec_in_ ./ sec_norm_;
    est_edc_sec_ = est_edc_sec(:, :, plotSecIdx(idxSec));
    est_edc_sec_ = est_edc_sec_ ./ sec_norm_;
    synth_edc_sec_ = synth_edc_sec(:, :, plotSecIdx(idxSec));
    synth_edc_sec_ = synth_edc_sec_ ./ sec_norm_;
    meas_edc_sec_out_ = meas_edc_sec_out(:, :, plotSecIdx(idxSec));
    meas_edc_sec_out_ = meas_edc_sec_out_ ./ sec_norm_;
    
    plot(time_axis, pow2db(meas_edc_sec_in_(:, plotBandIdx)), '-', 'linewidth', 2, 'SeriesIndex',idxSec);
    plot(time_axis, pow2db(est_edc_sec_(:, plotBandIdx)), ':', 'linewidth', 2, 'SeriesIndex',idxSec);
    %plot(time_axis, pow2db(synth_edc_sec_(:, plotBandIdx)), '-.', 'linewidth', 2, 'SeriesIndex',idxSec);
    plot(time_axis, pow2db(meas_edc_sec_out_(:, plotBandIdx)), '--', 'linewidth', 2, 'SeriesIndex',idxSec);
end
pc = colororder(gca());
for idxSec = 1:length(plotSecIdx)
    text(time_axis(end)-0.5+0.055*idxSec,-25, num2str(idxSec), 'Color',pc(idxSec, :));
end
grid on
ylim([-70, 0])
xlim([0, time_axis(end)])
xlabel("Time (s)")
ylabel("Magnitude (dB)")
legend("Measured", "Estimated", "Denoised")
%title("f="+ fbandsString(plotBandIdx))
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/EDC_fit"+num2str(target_snr)+".eps", 'Resolution',800)

%%% SECTOR EDC
if false  % Careful, could be many plots
% Sanity check: fitted EDC should match measured EDC
for secIdx=1:numSecs
figure;
hold on;
cmap = parula(numBands);
for bandIdx=1:numBands
    plot(time_axis, pow2db(meas_edc_sec(:, bandIdx, secIdx)), 'Color', cmap(bandIdx, :), 'LineWidth', 2, 'LineStyle', '-','DisplayName',fbandsString(bandIdx));
    plot(time_axis, pow2db(est_edc_sec(:, bandIdx, secIdx)), 'Color', cmap(bandIdx, :), 'LineWidth', 2, 'LineStyle', '--','DisplayName', "");
    plot(time_axis, pow2db(synth_edc_sec(:, bandIdx, secIdx)), 'Color', cmap(bandIdx, :), 'LineWidth', 2, 'LineStyle', ':','DisplayName', "");
end
xlabel('Time in s');
ylabel('Energy in dB');
legend()
%set(gca, 'FontSize', 15);
ylim([pow2db(max(meas_norm_edc_))-90, pow2db(max(meas_norm_edc_))+10]);
title("EDC meas vs. est (dashed) vs. synthesized (dotted), Sector " + num2str(secIdx))
end
end


%%% EDC Diff
% figure;
% sgtitle("EDC in dB at f="+ fbandsString(plotBandIdx))
% subplot(2, 1, 1)
% imagesc(pow2db(squeeze(meas_edc_sec(:, plotBandIdx, :))))
% subtitle("Measured EDC")
% %caxis([-100, 0])
% colorbar()
% 
% subplot(2, 1, 2)
% imagesc(pow2db(squeeze(est_edc_sec(:, plotBandIdx, :))))
% subtitle("Estimated EDC")
% xlabel("Sector Idx")
% xlabel("time in smp")
% ylabel("in dB")
% %caxis([-100, 0])
% colorbar()

% figure;
% sgtitle("Measured to Estimated EDC; Difference in dB")
% for bandIdx = 1:numBands
%     subplot(1, numBands+1, bandIdx)
%     imagesc(1:numSecs, time_axis, ...
%             10*log10(abs(squeeze(meas_edc_sec(1:fs, plotBandIdx, :)) - ...
%                          squeeze(est_edc_sec(1:fs, plotBandIdx, :)))));
%     set(gca,'YDir','normal')
%     subtitle("f="+ fbandsString(bandIdx))
%     caxis([-120, -60])
%     if bandIdx == 1; xlabel("Sector Idx"); ylabel("time in s"); end
% end
% ax = subplot(1, numBands+1,  numBands+1);
% axis off
% caxis([-120, -60])
% c = colorbar(ax);
% c.Label.String = 'dB';

% figure;
% sgtitle("Measured (solid) to Estimated (dashed) EDC")
% for bandIdx = 1:numBands
%     subplot(numBands, 1, bandIdx)
%     set(gca,'ColorOrderIndex',1)
%     plot(time_axis, ...
%          10*log10(abs(squeeze(meas_edc_sec(:, bandIdx, :)))), '-')
%     hold on
%     set(gca,'ColorOrderIndex',1)
%     plot(time_axis, ...
%         10*log10(abs(squeeze(est_edc_sec(:, bandIdx, :)))), '--')    
%     subtitle("f="+ fbandsString(bandIdx))
%     if bandIdx == numBands; xlabel("time in s"); ylabel("dB"); end
% end


%%% RMS
figure;
subplot(2, 1, 1)
stem(db(rms(rir_secs)))
hold on
stem(db(rms(rir_secs_synth)))
legend("Input", "Synthesized")
xlabel("Sector Idx")
ylabel("RMS in dB")
subplot(2, 1, 2)
stem(db(rms(rir_secs(lateIdx:end, :))))
hold on
stem(db(rms(rir_secs_synth(lateIdx:end, :))))
legend("Input Tail", "Synthesized Tail")
xlabel("Sector Idx")
ylabel("RMS in dB")

%%% Spec
[h_in,f] = freqz(sqrt(4*pi)*rir_secs_in_(lateIdx:end,1),1,2^12,fs);
[h_noisy,~] = freqz(sqrt(4*pi)*rir_secs(lateIdx:end,1),1,2^12,fs);
[h_synth,~] = freqz(sqrt(4*pi)*rir_secs_synth(lateIdx:end,1),1,2^12,fs);
[h_denoised,~] = freqz(sqrt(4*pi)*rir_secs_out_(lateIdx:end,1),1,2^12,fs);

figure; hold on; grid on;
%plot(w,h1);
%plot(w,h2);
plot(f,fracBandSmooth(mag2db(abs(h_in)), 3));
plot(f,fracBandSmooth(mag2db(abs(h_noisy)), 3));
plot(f,fracBandSmooth(mag2db(abs(h_synth)), 3));
plot(f,fracBandSmooth(mag2db(abs(h_denoised)), 3));
set(gca,'XScale','log');
xlim([20,fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Original Tail', 'Noisy Tail', 'Synth Tail', 'Denoised Tail');

figure;
hold on
grid on
for idxSec=1:length(plotSecIdx)
    plot3(zeros(size(f))+idxSec,f,...
            fracBandSmooth(mag2db(abs(h_secs_in(:,plotSecIdx(idxSec)))), 3),...
            '-', 'linewidth', 2, 'SeriesIndex',idxSec);
    plot3(zeros(size(f))+idxSec,f,...
            fracBandSmooth(mag2db(abs(h_secs_denoised(:,plotSecIdx(idxSec)))), 3),...
            ':', 'linewidth', 2, 'SeriesIndex',idxSec);
end
view(110,30)
set(gca, 'XDir','reverse')
set(gca,'YScale','log');
ylim([50,fs/2]);
xticks(1:1:length(plotSecIdx))
xlabel("idx")
ylabel("Frequency (Hz)")
zlabel("Magnitude (dB)")
%title("Tail Spectrum Noise-Free Input (solid) vs. Denoised (dashed)")
legend("Input", "Output", 'Location','northeast')
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/spec_secs_"+num2str(target_snr)+".eps", 'Resolution',800)

%%% Waterfall
figure;
[splot,fplot,tplot] = spectrogram(sqrt(4*pi)*rir_secs(lateIdx:end,1),...
                                  hanning(256),[],1024,fs);
meshz(tplot, fplot, mag2db(abs(2*splot)), 'facecolor', 'flat')  % 2: hann
caxis([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
view(110,30)
hold on
plot3(zeros(size(f))+0.001,f,...
        fracBandSmooth(mag2db(abs(h_noisy)), 3),...
        'k-', 'linewidth', 2);
title("Input Tail")
set(gca,'YScale','log');
ylim([50,fs/2]);
%zlim([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
zlim([-100, 20])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
zlabel("Magnitude (dB)")
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/water_in"+num2str(target_snr)+".eps", 'Resolution',800)

figure;
[splot,fplot,tplot] = spectrogram(sqrt(4*pi)*rir_secs_out_(lateIdx:end,1),...
                                      hanning(256),[],1024,fs);
meshz(tplot, fplot, mag2db(abs(2*splot)), 'facecolor', 'flat')  % 2: hann
caxis([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
view(110,30)
hold on
plot3(zeros(size(f))+0.001,f,...
        fracBandSmooth(mag2db(abs(h_denoised)), 3),...
        'k-', 'linewidth', 2);
title("Denoised Tail")
set(gca,'YScale','log');
ylim([50,fs/2]);
%zlim([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
zlim([-100, 20])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
zlabel("Magnitude (dB)")
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/water_out"+num2str(target_snr)+".eps", 'Resolution',800)


%%% SPH RMS
figure;
plotSphRms(rir_in_nm, 100, 0);
secPos = unitSph2cart(deg2rad(secDirs));
for idxSec = 1:size(secPos, 1)
          scatter(secDirs(idxSec, 1), secDirs(idxSec, 2), 30, 'black', 'filled',...
           'LineWidth', 1)
    if (any(idxSec == plotSecIdx))
        text(secDirs(idxSec, 1), secDirs(idxSec, 2),...
             num2str(idxSec), 'Color', 'white');
    else
         text(secDirs(idxSec, 1), secDirs(idxSec, 2),...
             num2str(idxSec), 'Color', '#838383');
    end
end
%title("Input RIR RMS")
%set(colorbar,'visible','off')
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/balloon.eps", 'Resolution',800)

% plotSphRms(rir_denoised_nm, 100);
% title("Denoised RIR RMS")
% 
% plotSphRms(rir_in_nm(lateIdx:end,:), 10000);
% title("Input RIR Tail RMS")
% plotSphRms(rir_noisy_nm(lateIdx:end,:), 10000);
% title("Input Noisy RIR Tail RMS")
% plotSphRms(rir_denoised_nm(lateIdx:end,:), 10000);
% title("Denoised RIR Tail RMS")

%%% SPH T60
figure;
[pltx, plty, pltz] = sph2cart(deg2rad(secDirs(:,1)), deg2rad(secDirs(:,2)), mean(measRT60sec,2));
scatter3(pltx, plty, pltz, 60, mean(measRT60sec,2), 'filled')
hold on
pltTri = delaunay(pltx, plty, pltz);
trisurf(pltTri,pltx, plty, pltz, mean(measRT60sec,2), "FaceAlpha",0.1, 'EdgeColor', 'none')
c=colorbar();
c.Label.String = "T_{60} (s)";
grid on
hold on
% axes
line([-1, 1], [0, 0], [0, 0], 'Color', 'black', 'LineWidth', 1)
line([0, 0], [-1, 1], [0, 0], 'Color', 'black', 'LineWidth', 1)
line([0, 0], [0, 0], [-1, 1], 'Color', 'black', 'LineWidth', 1)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis vis3d
view(-60, 25)


figure;
scatter(secDirs(:,1), secDirs(:,2), 60, mean(measRT60sec,2), 'filled')
c=colorbar();
c.Label.String = "T_{60} (s)";
set(gca, 'XDir','reverse')
xticks([-180:90:180])
xlim([-180, 180])
yticks([-90:45:90])
ylim([-90, 90])
xlabel('Azimuth')
ylabel('Elevation')
grid on
daspect([1 1 1])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/T60.jpg", 'Resolution',800)


figure
hold on
for idxSec = 1:length(plotSecIdx)
    plot([1:numBands], measRT60secIn(plotSecIdx(idxSec), :) + idxSec, 'o', 'LineWidth', 1.5, 'SeriesIndex',idxSec);
    plot([1:numBands], measRT60secOut(plotSecIdx(idxSec), :)+ idxSec, 'x', 'LineWidth', 1.5, 'SeriesIndex',idxSec);
    yline(idxSec, '--')
end
xticklabels(fbandsString)
legend("Input","Output", 'location', 'north')
xlabel("Frequency Band (Hz)")
ylabel("T_{60} + idx (s)")
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fig = gcf();set(fig, 'Units', 'Inches', 'Position', [0, 0, defFigsize(1), defFigsize(2)]);
exportgraphics(fig, "figs/compareT60_"+num2str(target_snr)+".eps", 'Resolution',800)

end %% PLOT


%% binauralize, listen to tail
if BIN
    rir_in_nm_bin = SHtoBin(rir_in_nm,fs);
    rir_synth_nm_bin = SHtoBin(rir_synth_nm_,fs);
    rir_denoised_nm_bin = SHtoBin(rir_denoised_nm,fs);
    rir_noisy_nm_bin = SHtoBin(rir_noisy_nm,fs);

    sound(5*rir_in_nm_bin(lateIdx:end,:), fs)
    pause(3)
    sound(5*rir_synth_nm_bin(lateIdx:end,:), fs)
    pause(3)
    sound(5*rir_denoised_nm_bin(lateIdx:end,:), fs)
    pause(3)
    sound(5*rir_noisy_nm_bin(lateIdx:end,:), fs)
else
    %careful, listening only to 0th order is only representative in full diffuse
    sound(10*sqrt(4*pi)*rir_in_nm(lateIdx:end,1), fs)
    pause(3)
    sound(10*sqrt(4*pi)*rir_synth_nm_(lateIdx:end,1), fs)
    pause(3)
    sound(10*sqrt(4*pi)*rir_denoised_nm(lateIdx:end,1), fs)
    pause(3)
    sound(10*sqrt(4*pi)*rir_noisy_nm(lateIdx:end,1), fs)
    pause(3)
end

%% Save
audiowrite("audio/RIRin"+target_snr+"N3D.wav", 0.8*rir_in_nm, fs,'BitsPerSample',32)
audiowrite("audio/RIRsynth"+target_snr+"N3D.wav", 0.8*rir_synth_nm_, fs,'BitsPerSample',32)
audiowrite("audio/RIRdenoised"+target_snr+"N3D.wav", 0.8*rir_denoised_nm, fs,'BitsPerSample',32)
audiowrite("audio/RIRnoisy"+target_snr+"N3D.wav", 0.8*rir_noisy_nm, fs,'BitsPerSample',32)


if BIN
    audiowrite("audio/RIRin"+target_snr+"Bin.wav", 1/4*rir_in_nm_bin, fs,'BitsPerSample',32)
    audiowrite("audio/RIRsynth"+target_snr+"Bin.wav", 1/4*rir_synth_nm_bin, fs,'BitsPerSample',32)
    audiowrite("audio/RIRdenoised"+target_snr+"Bin.wav", 1/4*rir_denoised_nm_bin, fs,'BitsPerSample',32)
    audiowrite("audio/RIRnoisy"+target_snr+"Bin.wav", 1/4*rir_noisy_nm_bin, fs,'BitsPerSample',32)

    % Render
    in_sig = audioread('audio/drums.wav');
    in_sig = cat(1, in_sig, zeros(numSmps, 1));
    rir_in_faded_bin = rir_in_nm_bin;
    %rir_in_faded_bin(end-1.*fs+1:end, :) = rir_in_faded_bin(end-1.*fs+1:end, :) .* db2mag(linspace(0, -30, 1*1.*fs)).';
    out_clean = fftfilt(rir_in_faded_bin, in_sig);
    out_noisy = fftfilt(rir_noisy_nm_bin, in_sig);
    out_denoised = fftfilt(rir_denoised_nm_bin, in_sig);
    audiowrite("audio/"+num2str(target_snr)+"clean.wav", out_clean/8, fs,'BitsPerSample',32)
    audiowrite("audio/"+num2str(target_snr)+"noisy.wav", out_noisy/8, fs,'BitsPerSample',32)
    audiowrite("audio/"+num2str(target_snr)+"denoised.wav", out_denoised/8, fs,'BitsPerSample',32)

end    



