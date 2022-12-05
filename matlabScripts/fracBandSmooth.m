function [h_smooth] = fracBandSmooth(h, frac)
% FRACBANDSMOOTH Pass magnitude h in dB, frac=1 for octave, frac=3 for 1/3.
% Chris Hold, 2022
    assert(isvector(h))
    h_smooth = zeros(size(h));
    for idx = 1:length(h)
        h_smooth(idx) = mean(h(max(floor(idx / 2^(1/(2*frac))), 1) : ...
                               min(ceil(idx * 2^(1/(2*frac))),length(h))));
    end
end
