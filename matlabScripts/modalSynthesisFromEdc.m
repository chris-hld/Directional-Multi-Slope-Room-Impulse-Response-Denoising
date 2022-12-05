function bandpassNoises = modalSynthesisFromEdc(edc,bands,numModes,fs)
len = length(edc);

thresh = 10^-30;
edc(edc < thresh) = thresh;
envelopes = (flipud(diff(flipud((edc)))));
envelopes = envelopes([1:end end],:);
envelopes(envelopes < thresh) = thresh;
envelopes = sqrt(envelopes);

% bands are the center frequencies, add bass and treble band
numberOfBands = length(bands)+2;

time = (1:len)';

% logspaced frequencies with a bit of jitter
modeFreq = logspace(log10(1),log10(fs/2),numModes+1);
modeFreq = modeFreq(1:end-1); % remove last because of randomizaton;

jitter = mean(diff(log10(modeFreq)));
modeFreq = 10.^(log10(modeFreq) + rand(1,numModes)*jitter);

modeFreq = modeFreq(modeFreq < fs/2 - 0.01*fs);  % Pull frequencies too close to aliasing
numModes = length(modeFreq);

% log space band index
bandFraction = (interp1(log10([1, bands, fs/2]),1:numberOfBands, log10(modeFreq)));
bandIndex = round(bandFraction);

bandFloor = floor(bandFraction);
bandMix = bandFraction - bandFloor;
% bandMix = round(bandMix);

modePhase = 2*pi*rand(1,numModes);

modalWeighting = sqrt(modeFreq/fs); % this is the stretching factor dg\df, where g is the logarithmic warping of the frequency (for example read https://siboehm.com/articles/19/normalizing-flow-network)

bandpassNoises = zeros(len,numberOfBands);

envelopesDB = mag2db(envelopes);

bandwiseModalWeighting = zeros(1,numberOfBands);
for it = 1:numModes
    % energy of sin is 1/sqrt(2)
    mode = sqrt(2) * sin(2*pi*time*modeFreq(it)/fs + modePhase(it));
    
    % modal weights which counteract the logarithmic spacing
    mode = mode .* modalWeighting(it);
    
    % apply edc envelopes
    mixEnvelope = db2mag((1-bandMix(it)) * (envelopesDB(:,bandFloor(it))) +...
                            bandMix(it) * (envelopesDB(:,bandFloor(it)+1)));
    mode = mode .* ( mixEnvelope );
    
    % sum to respectuve band
    bandpassNoises(:,bandIndex(it)) = bandpassNoises(:,bandIndex(it)) + mode;
end

for bit = 1:numberOfBands
    index = bandIndex == bit;
    bandwiseModalWeighting(bit) = sqrt(sum(modalWeighting(index).^2));
end


switch 2
    case 1
        % normalize the number of modes
        bandpassNoises = bandpassNoises / sqrt(numModes);

        % normalize the modal weighting broadband
        % sum of all bands has rms of 1 (without envelope)
        bandpassNoises = bandpassNoises / rms(modalWeighting);
    case 2  
        % normalize the modal weighting per band
        % each band has rms of 1 (without envelope)
        bandpassNoises = bandpassNoises ./ bandwiseModalWeighting;
    otherwise
        warning('not defined');
end

