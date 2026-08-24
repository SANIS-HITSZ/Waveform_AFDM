% Out-of-Band (OOB) emission
clear; clc; close all;
% rng(1);

color_set = [
    0.600, 0.400, 0.200;   % brown
    0.459, 0.439, 0.702;   % violet
    0.106, 0.620, 0.467;   % dark blue-green
    0.000, 0.447, 0.741;   % blue
    0.850, 0.325, 0.098    % red
];

%% Basic configurations
numRe        = 4096;                  % number of resource elements in resource grid
numRe_Rehd   = 4096;                  % Reheduled number of resource elements
lenCp        = 288;                   % length of prefix
rollOffWin   = 0.2;
numSym       = 200;                  % number of AFDM symbols
indxRe_Rehd  = [-numRe_Rehd/2: numRe_Rehd/2-1];     % Reheduled indices of resource elements
upSampCoef_time = 4;                % upsampling factor in the time domain
upSampCoef_freq = 8;                % spectrum upsampling
modulatOrder = 10;    % Reheduled modulation order

k_max = 3;
k_reserve = 4;     % reserved guard spacing for fractional doppler
chirpRate = (2*(ceil(k_max)+k_reserve)+1)/(2*numRe);  % chirp-rate
prechirpRate = 0;   % prechirp-rate

bitsRef = zeros(2^modulatOrder, modulatOrder);    % transmit bit-alphabet
qamsRef = zeros(2^modulatOrder, 1);               % transmit qam-alphabet
for indxQam = 0: 2^modulatOrder-1
    bitsRef(indxQam+1, :) = de2bi(indxQam, modulatOrder, 'left-msb');
    qamsRef(indxQam+1) = func_nrQamMapper(bitsRef(indxQam+1, :));
end

lenSegment = 2*numRe*upSampCoef_time;
lenOverlap = lenSegment / 2;
lenPSD     = numRe*upSampCoef_freq;
winAnalysis = hamming(lenSegment);

typeInterp = 1;
filterOrder  = numRe*upSampCoef_time;
rollOffTime  = 0.2; 
filterCoefTx = rootRaisedCosFilter(rollOffTime, filterOrder, upSampCoef_time);
filterCoefTx = filterCoefTx * sqrt(upSampCoef_time);

%% QAM data
symData = round(rand(numSym+2,numRe_Rehd)*2^modulatOrder - 0.5);   % transmit decimal symbols
qamData = qamsRef(symData+1);     % transmit qams
X = zeros(numSym+2, numRe);
X(:, mod(indxRe_Rehd,numRe)+1) = qamData;

%% CPP-AFDM
lenCp_ext = 0;
winTx = rectwin(numRe).';
s_overSamp = [];
for i_sym = 0:(numSym+2)-1
    % baseband
    s_sym = afdmModulator(X(i_sym+1,:), ...
        numRe, lenCp, lenCp_ext, winTx, chirpRate, prechirpRate);
    % upsampling
    if typeInterp == 0
        S_base = fftshift(fft(s_sym));
        S_interp = [zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2), ...
            S_base, zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2)];
        s_interp = ifft(ifftshift(S_interp)) * upSampCoef_time;
    else
        s_interp = zeros(1, upSampCoef_time*length(s_sym));
        s_interp(1: upSampCoef_time: end) = s_sym;
        s_interp = conv(s_interp, filterCoefTx);
        s_interp = s_interp(filterOrder/2+1: end-filterOrder/2);
    end
    %
    s_overSamp = [s_overSamp, s_interp];
end
s_overSamp = s_overSamp(upSampCoef_time*(numRe+lenCp+lenCp_ext)+1: ...
    end-upSampCoef_time*(numRe+lenCp+lenCp_ext));
% spectral estimation
[psd_cppAfdm, ~] = pwelch(s_overSamp, winAnalysis, lenOverlap, lenPSD, 'twosided');

%% PS-AFDM (legacy)
lenCp_ext = 0;
winTx = chebwin(numRe, 90).';
s_overSamp = [];
for i_sym = 0:(numSym+2)-1
    % baseband
    s_sym = afdmModulator(X(i_sym+1,:), ...
        numRe, lenCp, lenCp_ext, winTx, chirpRate, prechirpRate);
    % upsampling
    if typeInterp == 0
        S_base = fftshift(fft(s_sym));
        S_interp = [zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2), ...
            S_base, zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2)];
        s_interp = ifft(ifftshift(S_interp)) * upSampCoef_time;
    else
        s_interp = zeros(1, upSampCoef_time*length(s_sym));
        s_interp(1: upSampCoef_time: end) = s_sym;
        s_interp = conv(s_interp, filterCoefTx);
        s_interp = s_interp(filterOrder/2+1: end-filterOrder/2);
    end
    %
    s_overSamp = [s_overSamp, s_interp];
end
s_overSamp = s_overSamp(upSampCoef_time*(numRe+lenCp+lenCp_ext)+1: ...
    end-upSampCoef_time*(numRe+lenCp+lenCp_ext));
% spectral estimation
[psd_psAfdm, ~] = pwelch(s_overSamp, winAnalysis, lenOverlap, lenPSD, 'twosided');

%% WOLA/CAWOLA-AFDM
lenWin_oob = 8;
lenCp_ext = ceil(rollOffWin*numRe);   % length of extended cyclic-prefix for receiving windowing
lenCp_ext = lenCp_ext + mod(lenCp_ext,2);   % even number
winTx = rectwin(numRe).';
raiseEdgeLeft = 0.5*(1 - cos(pi*([0:lenWin_oob-1]+0.5)/lenWin_oob));
raiseEdgeRight = fliplr(raiseEdgeLeft);
winTx_oob = [raiseEdgeLeft, ones(1,numRe+lenCp+lenCp_ext-lenWin_oob), raiseEdgeRight];
s_overSamp = [];
for i_sym = 0:(numSym+2)-1
    % baseband
    s_sym = afdmModulator_oobSup(X(i_sym+1,:), ...
        numRe, lenCp, lenCp_ext, winTx, winTx_oob, chirpRate, prechirpRate);
    % upsampling
    if typeInterp == 0
        S_base = fftshift(fft(s_sym));
        S_interp = [zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2), ...
            S_base, zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2)];
        s_interp = ifft(ifftshift(S_interp)) * upSampCoef_time;
    else
        s_interp = zeros(1, upSampCoef_time*length(s_sym));
        s_interp(1: upSampCoef_time: end) = s_sym;
        s_interp = conv(s_interp, filterCoefTx);
        s_interp = s_interp(filterOrder/2+1: end-filterOrder/2);
    end
    %
    if i_sym == 0
        s_overSamp = [s_overSamp, s_interp];
    else
        s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
            = s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
                + s_interp(1: upSampCoef_time*lenWin_oob);
        s_overSamp = [s_overSamp, s_interp(upSampCoef_time*lenWin_oob+1: end)];
    end
end
s_overSamp = s_overSamp(upSampCoef_time*(numRe+lenCp+lenCp_ext)+1: ...
    end-upSampCoef_time*(numRe+lenCp+lenCp_ext+lenWin_oob));
% spectral estimation
[psd_wolaAfdm1, ~] = pwelch(s_overSamp, winAnalysis, lenOverlap, lenPSD, 'twosided');

%% WOLA/CAWOLA-AFDM
lenWin_oob = 16;
lenCp_ext = ceil(rollOffWin*numRe);   % length of extended cyclic-prefix for receiving windowing
lenCp_ext = lenCp_ext + mod(lenCp_ext,2);   % even number
winTx = rectwin(numRe).';
raiseEdgeLeft = 0.5*(1 - cos(pi*([0:lenWin_oob-1]+0.5)/lenWin_oob));
raiseEdgeRight = fliplr(raiseEdgeLeft);
winTx_oob = [raiseEdgeLeft, ones(1,numRe+lenCp+lenCp_ext-lenWin_oob), raiseEdgeRight];
s_overSamp = [];
for i_sym = 0:(numSym+2)-1
    % baseband
    s_sym = afdmModulator_oobSup(X(i_sym+1,:), ...
        numRe, lenCp, lenCp_ext, winTx, winTx_oob, chirpRate, prechirpRate);
    % upsampling
    if typeInterp == 0
        S_base = fftshift(fft(s_sym));
        S_interp = [zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2), ...
            S_base, zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2)];
        s_interp = ifft(ifftshift(S_interp)) * upSampCoef_time;
    else
        s_interp = zeros(1, upSampCoef_time*length(s_sym));
        s_interp(1: upSampCoef_time: end) = s_sym;
        s_interp = conv(s_interp, filterCoefTx);
        s_interp = s_interp(filterOrder/2+1: end-filterOrder/2);
    end
    %
    if i_sym == 0
        s_overSamp = [s_overSamp, s_interp];
    else
        s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
            = s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
                + s_interp(1: upSampCoef_time*lenWin_oob);
        s_overSamp = [s_overSamp, s_interp(upSampCoef_time*lenWin_oob+1: end)];
    end
end
s_overSamp = s_overSamp(upSampCoef_time*(numRe+lenCp+lenCp_ext)+1: ...
    end-upSampCoef_time*(numRe+lenCp+lenCp_ext+lenWin_oob));
% spectral estimation
[psd_wolaAfdm2, ~] = pwelch(s_overSamp, winAnalysis, lenOverlap, lenPSD, 'twosided');

%% WOLA/CAWOLA-AFDM
lenWin_oob = 32;
lenCp_ext = ceil(rollOffWin*numRe);   % length of extended cyclic-prefix for receiving windowing
lenCp_ext = lenCp_ext + mod(lenCp_ext,2);   % even number
winTx = rectwin(numRe).';
raiseEdgeLeft = 0.5*(1 - cos(pi*([0:lenWin_oob-1]+0.5)/lenWin_oob));
raiseEdgeRight = fliplr(raiseEdgeLeft);
winTx_oob = [raiseEdgeLeft, ones(1,numRe+lenCp+lenCp_ext-lenWin_oob), raiseEdgeRight];
s_overSamp = [];
for i_sym = 0:(numSym+2)-1
    % baseband
    s_sym = afdmModulator_oobSup(X(i_sym+1,:), ...
        numRe, lenCp, lenCp_ext, winTx, winTx_oob, chirpRate, prechirpRate);
    % upsampling
    if typeInterp == 0
        S_base = fftshift(fft(s_sym));
        S_interp = [zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2), ...
            S_base, zeros(1,(upSampCoef_time-1)*(numRe+lenCp+lenCp_ext)/2)];
        s_interp = ifft(ifftshift(S_interp)) * upSampCoef_time;
    else
        s_interp = zeros(1, upSampCoef_time*length(s_sym));
        s_interp(1: upSampCoef_time: end) = s_sym;
        s_interp = conv(s_interp, filterCoefTx);
        s_interp = s_interp(filterOrder/2+1: end-filterOrder/2);
    end
    %
    if i_sym == 0
        s_overSamp = [s_overSamp, s_interp];
    else
        s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
            = s_overSamp(end-upSampCoef_time*lenWin_oob+1: end) ...
                + s_interp(1: upSampCoef_time*lenWin_oob);
        s_overSamp = [s_overSamp, s_interp(upSampCoef_time*lenWin_oob+1: end)];
    end
end
s_overSamp = s_overSamp(upSampCoef_time*(numRe+lenCp+lenCp_ext)+1: ...
    end-upSampCoef_time*(numRe+lenCp+lenCp_ext+lenWin_oob));
% spectral estimation
[psd_wolaAfdm3, ~] = pwelch(s_overSamp, winAnalysis, lenOverlap, lenPSD, 'twosided');

%% Figures
powerOffset_welch = 2*pi/upSampCoef_time;
freqLabel = [-lenPSD/2: lenPSD/2-1] * (upSampCoef_time*numRe)/lenPSD;
figure; 
plot(freqLabel, 10*log10(fftshift(psd_cppAfdm*powerOffset_welch)),'Color',color_set(1,:),'linewidth',1);
hold on;
plot(freqLabel, 10*log10(fftshift(psd_psAfdm*powerOffset_welch)),'Color',color_set(2,:),'linewidth',1);
plot(freqLabel, 10*log10(fftshift(psd_wolaAfdm1*powerOffset_welch)),'Color',color_set(3,:),'linewidth',1);
plot(freqLabel, 10*log10(fftshift(psd_wolaAfdm2*powerOffset_welch)),'Color',color_set(4,:),'linewidth',1);
plot(freqLabel, 10*log10(fftshift(psd_wolaAfdm3*powerOffset_welch)),'Color',color_set(5,:),'linewidth',1);
xline(-numRe*(1+rollOffTime)/2, 'k:', 'LineWidth', 1.5);
xline(numRe*(1+rollOffTime)/2-1, 'k:', 'LineWidth', 1.5);
xlabel('Normalized Frequency $f/\Delta_{\mathrm{F}}$', 'interpreter', 'latex', 'FontSize', 12);
ylabel('Normalized Magnitude (dB)', 'interpreter', 'latex', 'FontSize', 12);
legend({'CPP-AFDM [7]', 'PS-AFDM [10] (Tx. Cheb.90dB Win.)', ...
    '(proposed) (CA)WOLA-AFDM (Tx. RC. Win. $L_{\mathrm{O}}=8$)', ...
    '(proposed) (CA)WOLA-AFDM (Tx. RC. Win. $L_{\mathrm{O}}=16$)', ...
    '(proposed) (CA)WOLA-AFDM (Tx. RC. Win. $L_{\mathrm{O}}=32$)'}, ...
    'interpreter', 'latex', 'Box', 'off', 'FontSize', 10);
grid on;
axis([0, 2*numRe, -90, 10]);


%% AFDM modulator
function s_afdm = afdmModulator(xData, ...
    numRe, lenCp, lenCp_ext, winTx, chirpRate, prechirpRate)
    % chirp base function
    chirpBase1 = exp(1i * 2 * pi * chirpRate * [-(lenCp+lenCp_ext): numRe-1].^2);
    chirpBase2 = exp(1i * 2 * pi * prechirpRate * [0: numRe-1].^2);
    % prechirping and ifft
    s_ofdm = ifft(xData .* chirpBase2) * sqrt(numRe);
    % transmit shaping windowing
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    s_ofdmW = s_ofdm .* winTx;
    % cyclic-prefix attachment
    s_ofdmCp = [s_ofdmW(end-(lenCp+lenCp_ext)+1:end), s_ofdmW];
    % chirping
    s_afdm = s_ofdmCp .* chirpBase1;
end

function s_afdm = afdmModulator_oobSup(xData, ...
    numRe, lenCp, lenCp_ext, winTx, winTx_oob, chirpRate, prechirpRate)
    % chirp base function
    lenWin_oob = length(winTx_oob) - (numRe+lenCp+lenCp_ext);
    chirpBase1 = exp(1i * 2 * pi * chirpRate * [-(lenCp+lenCp_ext): numRe+lenWin_oob-1].^2);
    chirpBase2 = exp(1i * 2 * pi * prechirpRate * [0: numRe-1].^2);
    % prechirping and ifft
    s_ofdm = ifft(xData .* chirpBase2) * sqrt(numRe);
    % transmit shaping windowing
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    s_ofdmW = s_ofdm .* winTx;
    % cyclic-prefix attachment
    s_ofdmCp = [s_ofdmW(end-(lenCp+lenCp_ext)+1:end), s_ofdmW];
    % oob suppress windowing
    s_ofdmCp = [s_ofdmCp, s_ofdmW(1:lenWin_oob)];
    s_ofdmCp = s_ofdmCp .* winTx_oob;
    % chirping
    s_afdm = s_ofdmCp .* chirpBase1;
end

%% Root raised cosine pulse shaping filter
function filterCoef = rootRaisedCosFilter(rollOff, filterOrder, upSampCoef)
    % Coded by Haojian Zhang
    % UWB-LAB, Rehool of Information Reience and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    indxFreq = [-filterOrder/2: filterOrder/2] / (filterOrder+1);
    winFilter = zeros(1, filterOrder+1);
    winFilter(abs(indxFreq) <= (1-rollOff)/(2*upSampCoef)) = 1;
    indxFreq_other = indxFreq((abs(indxFreq)>(1-rollOff)/(2*upSampCoef)) ...
        &(abs(indxFreq)<=(1+rollOff)/(2*upSampCoef)));
    winFilter((abs(indxFreq)>(1-rollOff)/(2*upSampCoef)) ...
        &(abs(indxFreq)<=(1+rollOff)/(2*upSampCoef))) ...
        = cos(pi/(2*rollOff)*(abs(indxFreq_other)*upSampCoef-(1-rollOff)/2)).^2;
    winFilter = sqrt(winFilter);
    filterCoef = ifftshift(ifft(ifftshift(winFilter)));
    filterCoef = filterCoef / sqrt(sum(abs(filterCoef).^2));
end

%% QAM-mapper
function qam = func_nrQamMapper(bits)
    % Bit-to-QAM Mapping in 5G-NR (BPSK is not included here)
    % Ref: 3GPP 38.211
    %
    % Input:
    %     bits,  bits to be mapped to a qam
    % Output:
    %     qam,   a mapped qam
    %
    % Coded by Haojian Zhang
    % UWB-LAB, Rehool of Information Reience and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    modulatOrder = length(bits);  % Reheduled modulation order
    switch modulatOrder
        case 2  % qpsk
            qam =      (1 - 2*bits(1)) ...
                + 1i * (1 - 2*bits(2));
            qam = qam / sqrt(2);
        case 4  % 16qam
            qam =      (1 - 2*bits(1)) * (2 - (1 - 2*bits(3))) ...
                + 1i * (1 - 2*bits(2)) * (2 - (1 - 2*bits(4)));
            qam = qam / sqrt(10);
        case 6  % 64qam
            qam =      (1 - 2*bits(1)) * (4 - (1 - 2*bits(3)) * (2 - (1 - 2*bits(5)))) ...
                + 1i * (1 - 2*bits(2)) * (4 - (1 - 2*bits(4)) * (2 - (1 - 2*bits(6))));
            qam = qam / sqrt(42);
        case 8   % 256qam
            qam =      (1 - 2*bits(1)) * (8 - (1 - 2*bits(3)) * (4 - (1 - 2*bits(5)) * (2 - (1 - 2*bits(7))))) ...
                + 1i * (1 - 2*bits(2)) * (8 - (1 - 2*bits(4)) * (4 - (1 - 2*bits(6)) * (2 - (1 - 2*bits(8)))));
            qam = qam / sqrt(170);
        case 10   % 1024qam
            qam =      (1 - 2*bits(1)) * (16 - (1 - 2*bits(3)) * (8 - (1 - 2*bits(5)) * (4 - (1 - 2*bits(7)) * (2 - (1 - 2*bits(9)))))) ...
                + 1i * (1 - 2*bits(2)) * (16 - (1 - 2*bits(4)) * (8 - (1 - 2*bits(6)) * (4 - (1 - 2*bits(8)) * (2 - (1 - 2*bits(10))))));
            qam = qam / sqrt(682);
    end
end
