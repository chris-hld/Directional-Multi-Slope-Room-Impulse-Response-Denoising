function [SNR_diff,specDiffs,rt60diff] = compareMetrics(rir_in_nm,rir_out_nm,edcs,pars,PLOT)
%COMPAREMETRICS Calculate and plot some input to output metrics.
% Uses the output of DIRECTIONAL_DENOISE_SRIR.m
N_sph = sqrt(length(rir_in_nm(1,:)))-1; % SH order
numSmps = size(rir_in_nm, 1);
fs = pars.fs;
numBands = pars.numBands;
timeAxis = linspace(0, (numSmps - 1) / fs, numSmps );

measEdcSec = edcs(:,:,:,1);
estEdcSec = edcs(:,:,:,2);
synthEdcSec = edcs(:,:,:,3);


% Onset
[~, k_max] = max(abs(rir_in_nm(:, 1)));
k_max = k_max - 10;  % get the first couple of onset samples
lateIdx = k_max + round(0.1 * fs);

% Decompose again
[~, secDirs] = getTdesign(2*N_sph);
secPos = unitSph2cart(secDirs);
R = calculateRotationMatrix(secPos(1, :), [1, 0, 0]);
secPos = secPos * R;  % rotate to sec0 = [0,0]
[secAzi, secEle, ~] = cart2sph(secPos(:,1), secPos(:,2), secPos(:,3));
secDirs = [rad2deg(secAzi), rad2deg(secEle)];
numSecs = size(secDirs, 1);

[A, B_AP] = designSphFilterBank(N_sph,secDirs,pars.spatFilterCoeffs,'AP');
rirSecsIn_ = rir_in_nm * A;
rirSecsOut_ = rir_out_nm * A;

%% RT60
secRT60 = zeros(1, numSecs);
measRT60sec = zeros(numSecs, numBands);
estRT60sec = zeros(numSecs, numBands);
for idxSec = 1:numSecs
    meas_edc_sec_ = measEdcSec(:, :, idxSec);
    meas_edc_sec_ = meas_edc_sec_ ./ max(meas_edc_sec_, [], 1);
    est_edc_sec_ = estEdcSec(:, :, idxSec);
    est_edc_sec_ = est_edc_sec_ ./ max(est_edc_sec_, [], 1);
    for idxBand = 1:numBands
        measRT60sec(idxSec,idxBand) = 2*(find(meas_edc_sec_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
        estRT60sec(idxSec,idxBand) = 2*(find(est_edc_sec_(:,idxBand) > db2pow(-30), 1, 'last') ./ fs);
    end
end


measEdcSecIn = zeros(size(measEdcSec));
measEdcSecOut = zeros(size(measEdcSec));
measRT60secIn = zeros(numSecs, numBands);
measRT60secOut = zeros(numSecs, numBands);
for idxSec = 1:numSecs
    [meas_edc_sec_in_, meas_norm_edc_sec_in_] = rir2decay(rirSecsIn_(:,idxSec), fs, pars.fBands, true, true, true, pars.includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    [meas_edc_sec_out_, meas_norm_edc_sec_out_] = rir2decay(rirSecsOut_(:,idxSec), fs, pars.fBands, true, true, true, pars.includeResidualBands); % bools: doBackwardsInt, analyseFullRIR, normalize, includeResidualBands
    measEdcSecIn(:, :, idxSec) = meas_edc_sec_in_ .* meas_norm_edc_sec_in_;
    measEdcSecOut(:, :, idxSec) = meas_edc_sec_out_ .* meas_norm_edc_sec_out_;
    
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
rt60diff = measRT60secOut - measRT60secIn;

%% SNR
rangeSNR = round(0.1 * fs);
% 4 pi might be missing here, but should cancel
SNR_in = trace(rir_in_nm(k_max:k_max+rangeSNR, :)'*rir_in_nm(k_max:k_max+rangeSNR, :)) / ...
    trace(rir_in_nm(end-rangeSNR:end, :)'*rir_in_nm(end-rangeSNR:end, :));
SNR_out = trace(rir_out_nm(k_max:k_max+rangeSNR, :)'*rir_out_nm(k_max:k_max+rangeSNR, :)) / ...
    trace(rir_out_nm(end-rangeSNR:end, :)'*rir_out_nm(end-rangeSNR:end, :));
SNR_diff = 10*log10(SNR_out/SNR_in);

%% SPEC
% Range should be adjusted if not a clean SRIR to exclude noise floor
rangeSpec = round(1 * fs);

hSecsIn = zeros(4096,numSecs);
hSecsOut = zeros(4096,numSecs);
specDiffs = zeros(4096,numSecs);
for idxSec=1:numSecs
    [hSecsIn(:,idxSec),f] = freqz(rirSecsIn_(lateIdx:lateIdx+rangeSpec,idxSec),1,2^12,fs);
    [hSecsOut(:,idxSec),~] = freqz(rirSecsOut_(lateIdx:lateIdx+rangeSpec,idxSec),1,2^12,fs);
    specDiffs(:,idxSec) = fracBandSmooth(mag2db(abs(hSecsOut(:,idxSec))), 3) - ...
        fracBandSmooth(mag2db(abs(hSecsIn(:,idxSec))), 3);
end



if PLOT
    disp("Plotting...")
%defFigsize = [3.25, 2.5];
plotBandIdx = 5;
%plotSecIdx = [1, 8, 11, 15, 18, 23];
plotSecIdx = 1:4:numSecs;
if pars.includeResidualBands
    fbandsString = ["Lo" num2cell(pars.fBands) "Hi"];
else
    fbandsString = string(num2cell(pars.fBands));
end


% IR
figure;
hold on;
plot(timeAxis, pow2db(dot(rir_in_nm, rir_in_nm, 2)), ':');
plot(timeAxis, pow2db(dot(rir_out_nm, rir_out_nm, 2)), '-');
xlabel("Time (s)")
ylabel("Magnitude (dB)")
legend("Input", "Output")
grid on
title("Signal Energy")

% EDC
figure;
hold on
%sec_norm_ = max(meas_edc_sec(:, :, 1), [], 1);  % norm on 1st sec
for idxSec = 1:length(plotSecIdx)
    meas_edc_sec_in_ = measEdcSecIn(:, :, plotSecIdx(idxSec));
    sec_norm_ = max(meas_edc_sec_in_, [], 1);
    meas_edc_sec_in_ = meas_edc_sec_in_ ./ sec_norm_;
    est_edc_sec_ = estEdcSec(:, :, plotSecIdx(idxSec));
    est_edc_sec_ = est_edc_sec_ ./ sec_norm_;
    %synth_edc_sec_ = synthEdcSec(:, :, plotSecIdx(idxSec));
    %synth_edc_sec_ = synth_edc_sec_ ./ sec_norm_;
    meas_edc_sec_out_ = measEdcSecOut(:, :, plotSecIdx(idxSec));
    meas_edc_sec_out_ = meas_edc_sec_out_ ./ sec_norm_;
    
    plot(timeAxis, pow2db(meas_edc_sec_in_(:, plotBandIdx)), '-', 'linewidth', 2, 'SeriesIndex',idxSec);
    plot(timeAxis, pow2db(est_edc_sec_(:, plotBandIdx)), ':', 'linewidth', 2, 'SeriesIndex',idxSec);
    %plot(timeAxis, pow2db(synth_edc_sec_(:, plotBandIdx)), '-.', 'linewidth', 2, 'SeriesIndex',idxSec);
    plot(timeAxis, pow2db(meas_edc_sec_out_(:, plotBandIdx)), '--', 'linewidth', 2, 'SeriesIndex',idxSec);
end
pc = colororder(gca());
for idxSec = 1:length(plotSecIdx)
    text(timeAxis(end)-0.5+0.055*idxSec,-25, num2str(plotSecIdx(idxSec)), 'Color',pc(idxSec, :));
end
grid on
ylim([-70, 0])
xlim([0, timeAxis(end)])
xlabel("Time (s)")
ylabel("Magnitude (dB)")
legend("Measured", "Estimated", "Denoised")
title("Directional Energy Decay")

% Waterfall
figure;
[splot,fplot,tplot] = spectrogram(rirSecsIn_(lateIdx:end,1),...
                                  hanning(256),[],1024,fs);
meshz(tplot, fplot, mag2db(abs(2*splot)), 'facecolor', 'flat')  % 2: hann
view(110,30)
hold on
plot3(zeros(size(f))+0.001,f,...
        fracBandSmooth(mag2db(abs(hSecsIn(:,1))), 3),...
        'k-', 'linewidth', 2);
title("Input Tail")
set(gca,'YScale','log');
ylim([50,fs/2]);
%zlim([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
zlim([-100, 20])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
zlabel("Magnitude (dB)")

figure;
[splot,fplot,tplot] = spectrogram(rirSecsOut_(lateIdx:end,1),...
                                      hanning(256),[],1024,fs);
meshz(tplot, fplot, mag2db(abs(2*splot)), 'facecolor', 'flat')  % 2: hann
view(110,30)
hold on
plot3(zeros(size(f))+0.001,f,...
        fracBandSmooth(mag2db(abs(hSecsOut(:,1))), 3),...
        'k-', 'linewidth', 2);
title("Output Tail")
set(gca,'YScale','log');
ylim([50,fs/2]);
%zlim([max(mag2db(abs(h_in)))-120, max(mag2db(abs(h_in)))])
zlim([-100, 20])
xlabel("Time (s)")
ylabel("Frequency (Hz)")
zlabel("Magnitude (dB)")

% SPH RMS
plotSphRms(rir_in_nm, 1, 1, [], secDirs);

end %PLOT

end