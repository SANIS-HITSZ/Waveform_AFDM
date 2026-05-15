% Simulation of the proposed OS-PS-AFDM transceiver
% Coded by Haojian Zhang
% UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
% Copyright (c) 2025, all rights reserved.
clc; clear; close all;
rng(1013);

%% Parameters of signalling
numSc = 4096;         % number of subcarriers in resource grid
lenCp = 288;          % length of normal cyclic-prefix
modulatOrder = 10;    % scheduled modulation order
filterOrder = 8192;   % order of the time-domain pulse shaping filter (after up-sampling)
rollOffTime = 0.2;    % roll-off factor for the time-domain pulse shaping (root-raised-cosine)
upSampleCoef = 8;     % upsampling factor in the time domain
lenPulseTime = 32;    % considered valid length of the overall time-domain shaping pulse
bitsRef = zeros(2^modulatOrder, modulatOrder);    % transmit bit-alphabet
qamsRef = zeros(1, 2^modulatOrder);               % transmit qam-alphabet
for indxQam = 0: 2^modulatOrder-1
    bitsRef(indxQam+1, :) = de2bi(indxQam, modulatOrder, 'left-msb');
    qamsRef(indxQam+1) = func_nrQamMapper(bitsRef(indxQam+1, :));
end
startSc_schd = 0;     % scheduled start subcarrier index
numSc_schd = 50*12;   % scheduled subcarrier number
indxSc_schd = [startSc_schd: startSc_schd + numSc_schd - 1];  % scheduled subcarrier indexes

%% Monte Carlo tests
numMcTest = 10000;                   % the number of monte carlo tests
rollOffWin_set = [0: 0.05: 0.4];    % roll-off factor for receiving windowing
sinr_set = [10: 5: 35];             % received snr (db)

ber_imperfCSI_all = zeros(length(rollOffWin_set)+2, length(sinr_set), numMcTest);
ber_perfCSI_all = zeros(length(rollOffWin_set)+2, length(sinr_set), numMcTest);
condnum_all = zeros(length(rollOffWin_set)+2, numMcTest);
nmse_all = zeros(length(rollOffWin_set)+2, numMcTest);

for indxMcTest = 0: numMcTest-1

    %% Parameters of delay-Doppler channel
    numPath = 10;   % paht number
    k_max = 3;      % maximum doppler
    ell_max = 10;   % maximum random range of delay
    tapDelayPath = rand(1,numPath)*ell_max + lenPulseTime/2;  % delay taps
    tapDopplerPath = cos(rand(1,numPath)*2*pi)*k_max;         % doppler taps
    tapGainPath = (randn(1,numPath)+1i*randn(1,numPath))/sqrt(numPath*2);  % rayleigh fading path gains

    %% Chirp-rate
    k_reserve = 4;     % reserved guard spacing for fractional doppler
    chirpParam1 = (2*(k_max+k_reserve)+1)/(2*numSc);  % chirp-rate
    chirpParam2 = 0;   % prechirp-rate

    %% Data modulation symbols
    symData_user = round(rand(1, numSc_schd)*2^modulatOrder - 0.5);   % transmit decimal symbols
    bitData_user = bitsRef(symData_user+1, :);  % transmit bits
    xData_user   = qamsRef(symData_user+1);     % transmit qams

    %% Resource grid mapping
    xData = zeros(1, numSc);                 % resource grid
    xData(indxSc_schd+1) = xData_user;       % data mapping
    xData = xData * sqrt(numSc/numSc_schd);  % power boosting
    indxPilot = indxSc_schd(end);            % index of embedded pilot (indivisual pilot)
    xPilot = zeros(1,numSc);                 % resource grid
    xPilot(indxPilot+1) = 1;                 % pilot mapping

    %% Additive white noise
    lenCp_ext = ceil(max(rollOffWin_set)*numSc);   % maximum length of extended cyclic-prefix for receiving windowing
    lenCp_ext = lenCp_ext + mod(lenCp_ext,2);
    numSample = numSc + lenCp + lenCp_ext;
    a_wgn = (randn(1,numSample)+1i*randn(1,numSample)) / sqrt(2);  % normalized additive white noise

    %% Different shaping windows
    for indxWin = 0: (length(rollOffWin_set)+2)-1

        %% Shaping window
        if indxWin < length(rollOffWin_set)
            % proposed
            rollOffWin = rollOffWin_set(indxWin+1);
            lenCp_ext = ceil(rollOffWin*numSc);   % length of extended cyclic-prefix for receiving windowing
            lenCp_ext = lenCp_ext + mod(lenCp_ext,2);   % even number
            rollOffWin = lenCp_ext / numSc;   % roll-off factor modification
            winTx = rectwin(numSc).';   % transmit shaping window
            winRx = raisedCosWin(rollOffWin, numSc);   % receive shaping window
        elseif indxWin == length(rollOffWin_set)
            % legacy
            lenCp_ext = 0;   % length of extended cyclic-prefix for receiving windowing
            winTx = chebwin(numSc, 70).';   % transmit shaping window
            winRx = rectwin(numSc).';       % receive shaping window
        elseif indxWin == length(rollOffWin_set)+1
            % legacy
            lenCp_ext = 0;   % length of extended cyclic-prefix for receiving windowing
            winTx = chebwin(numSc, 90).';   % transmit shaping window
            winRx = rectwin(numSc).';       % receive shaping window
        end

        %% AFDM transceiver procedures
        numSample = numSc + lenCp + lenCp_ext;
        a_wgn_temp = a_wgn(1: numSample);
        % data
        [yData, yNoise] = afdmTransceiver(xData, a_wgn_temp, ...
            numSc, lenCp, lenCp_ext, winTx, winRx, ...
            chirpParam1, chirpParam2, ...
            rollOffTime, filterOrder, upSampleCoef, ...
            numPath, tapGainPath, tapDelayPath, tapDopplerPath);
        % pilot
        [yPilot, ~] = afdmTransceiver(xPilot, a_wgn_temp, ...
            numSc, lenCp, lenCp_ext, winTx, winRx, ...
            chirpParam1, chirpParam2, ...
            rollOffTime, filterOrder, upSampleCoef, ...
            numPath, tapGainPath, tapDelayPath, tapDopplerPath);

        %% Channel matrix construction in noiseless condition
        % perfect channel matrix
        if indxWin < length(rollOffWin_set)
            hData = inOutAfdm_rxWin(numSc, indxSc_schd, lenCp_ext, ...
                chirpParam1, chirpParam2, rollOffTime, lenPulseTime, winRx, ...
                numPath, tapGainPath, tapDelayPath, tapDopplerPath);
            % condition number
            [~, S, ~] = svd(hData'*hData);
            S = sqrt(abs(diag(S)));
            condnum_all(indxWin+1, indxMcTest+1) = max(S) / min(S);
        else
            hData = inOutAfdm_txWin(numSc, indxSc_schd, ...
                chirpParam1, chirpParam2, rollOffTime, lenPulseTime, winTx, ...
                numPath, tapGainPath, tapDelayPath, tapDopplerPath);
            % condition number
            [~, S, ~] = svd(hData'*hData);
            S = sqrt(abs(diag(S)));
            condnum_all(indxWin+1, indxMcTest+1) = max(S) / min(S);
        end
        % estimated channel matrix
        numDoppler = abs(2 * numSc * chirpParam1);
        numDelay = ceil(ell_max) + lenPulseTime + 1;
        labelDoppler = [-(numDoppler-1)/2: (numDoppler-1)/2];
        labelDelay = [0: ceil(ell_max)+lenPulseTime];
        hMat_est = chanMat_est(yPilot, xPilot, indxPilot, ...
            numSc, indxSc_schd, chirpParam1, chirpParam2, labelDelay, labelDoppler);
        % nmse under noiseless conditions
        yData_est = (hMat_est * (xData_user*sqrt(numSc/numSc_schd)).').';
        nmse_all(indxWin+1, indxMcTest+1) = sum(abs(yData-yData_est).^2) / sum(abs(yData).^2);

        %% Different SINRs
        for indxSinr = 0: length(sinr_set)-1

            %% Noisy received symbols
            powerNoise = 10^(-sinr_set(indxSinr+1)/10);
            yData_noise = (yData + sqrt(powerNoise)*yNoise).';

            %% Power de-boosting
            yData_noise = yData_noise / sqrt(numSc/numSc_schd);
            powerNoise = powerNoise / (numSc/numSc_schd);

            %% LMMSE equalization
            xData_perfCSI   = (hData'   *hData    + powerNoise*eye(numSc_schd)) \ (hData'   *yData_noise);
            xData_imperfCSI = (hMat_est'*hMat_est + powerNoise*eye(numSc_schd)) \ (hMat_est'*yData_noise);
            bitEst_perfCSI   = zeros(numSc_schd, modulatOrder);
            bitEst_imperfCSI = zeros(numSc_schd, modulatOrder);
            for indxX = 0: numSc_schd-1
                [~, indxMin] = min(abs(xData_perfCSI(indxX+1) - qamsRef));
                bitEst_perfCSI(indxX+1, :) = bitsRef(indxMin, :);
                [~, indxMin] = min(abs(xData_imperfCSI(indxX+1) - qamsRef));
                bitEst_imperfCSI(indxX+1, :) = bitsRef(indxMin, :);
            end
            ber_perfCSI_all(indxWin+1, indxSinr+1, indxMcTest+1) = ...
                sum(sum(bitEst_perfCSI   ~= bitData_user)) / (numSc_schd*modulatOrder);
            ber_imperfCSI_all(indxWin+1, indxSinr+1, indxMcTest+1) = ...
                sum(sum(bitEst_imperfCSI ~= bitData_user)) / (numSc_schd*modulatOrder);

        end
        
        disp([num2str(indxWin+1),' / ',num2str(length(rollOffWin_set)+2)]);

    end

    disp([num2str(indxMcTest+1),' / ',num2str(numMcTest)]);

end

save('data_ber_imperfCSI.mat', 'ber_imperfCSI_all');
save('data_ber_perfCSI.mat', 'ber_perfCSI_all');
save('data_condnum.mat', 'condnum_all');
save('data_NMSE.mat', 'nmse_all');


%% Figures
marker_set = {'+', 'o', '*', '.', 's', '>', '<', 'p', 'd'};
% BER using imperfect CSI
label_set = {};
figure;
for indxWin = 0: (length(rollOffWin_set)+2)-1
    if indxWin < length(rollOffWin_set)
        semilogy(sinr_set, mean(ber_imperfCSI_all(indxWin+1, :, :), 3), ['r-',marker_set{indxWin+1}], 'lineWidth', 1);
        hold on;
        rollOffWin = rollOffWin_set(indxWin+1);
        label_str = {['(proposed) Rx. RC Win.\&Overlap, $\alpha_{\mathrm{W}}=$', num2str(rollOffWin)]};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)
        semilogy(sinr_set, mean(ber_imperfCSI_all(indxWin+1, :, :), 3), 'b--+', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 70dB-sidelobes'};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)+1
        semilogy(sinr_set, mean(ber_imperfCSI_all(indxWin+1, :, :), 3), 'b--o', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 90dB-sidelobes'};
        label_set = [label_set, label_str];
    end
end
xlabel('Received SNR (dB)', 'fontsize', 12, 'interpreter', 'latex');
ylabel('Uncoded BER (Estimated $\hat{\mathbf{H}}$, Noiseless)', 'fontsize', 12, 'interpreter', 'latex');
legend(label_set, 'fontsize', 10, 'interpreter', 'latex', 'Box', 'Off');
grid on;
ax = gca;
ax.TickLabelInterpreter = 'latex';
axis([min(sinr_set), max(sinr_set), 1e-4, 1]);
set(gcf, 'PaperType', 'a3');
% BER using perfect CSI
label_set = {};
figure;
for indxWin = 0: (length(rollOffWin_set)+2)-1
    if indxWin < length(rollOffWin_set)
        semilogy(sinr_set, mean(ber_perfCSI_all(indxWin+1, :, :), 3), ['r-',marker_set{indxWin+1}], 'lineWidth', 1);
        hold on;
        rollOffWin = rollOffWin_set(indxWin+1);
        label_str = {['(proposed) Rx. RC Win.\&Overlap, $\alpha_{\mathrm{W}}=$', num2str(rollOffWin)]};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)
        semilogy(sinr_set, mean(ber_perfCSI_all(indxWin+1, :, :), 3), 'b--+', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 70dB-sidelobes'};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)+1
        semilogy(sinr_set, mean(ber_perfCSI_all(indxWin+1, :, :), 3), 'b--o', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 90dB-sidelobes'};
        label_set = [label_set, label_str];
    end
end
xlabel('Received SNR (dB)', 'fontsize', 12, 'interpreter', 'latex');
ylabel('Uncoded BER (Perfect $\mathbf{H}$)', 'fontsize', 12, 'interpreter', 'latex');
legend(label_set, 'fontsize', 10, 'interpreter', 'latex', 'Box', 'Off');
grid on;
ax = gca;
ax.TickLabelInterpreter = 'latex';
axis([min(sinr_set), max(sinr_set), 1e-4, 1]);
set(gcf, 'PaperType', 'a3');
% nmse floor of channel estimation
figure;
plot(rollOffWin_set, 10*log10(mean(nmse_all(1: length(rollOffWin_set),:), 2)), 'r', 'lineWidth', 1);
hold on;
plot(rollOffWin_set, ones(1,length(rollOffWin_set))*10*log10(mean(nmse_all(length(rollOffWin_set)+1,:), 2)), 'b--+', 'lineWidth', 1);
plot(rollOffWin_set, ones(1,length(rollOffWin_set))*10*log10(mean(nmse_all(length(rollOffWin_set)+2,:), 2)), 'b--o', 'lineWidth', 1);
xlabel('Roll-Off Factor $\alpha_{\mathrm{W}}$', 'fontsize', 12, 'interpreter', 'latex');
ylabel('NMSE Floor of Estimated $\hat{\mathbf{H}}$ (dB)', 'fontsize', 12, 'interpreter', 'latex');
legend({'(proposed) Rx. RC Win.\&Overlap', ...
        '(legacy) Tx. Chebyshev Win., 70dB-sidelobes', ...
        '(legacy) Tx. Chebyshev Win., 90dB-sidelobes'}, ...
        'fontsize', 10, 'interpreter', 'latex', 'Box', 'Off');
grid on;
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gcf, 'PaperType', 'a3');
% condition number of channel matrix
label_set = {};
figure;
for indxWin = 0: (length(rollOffWin_set)+2)-1
    data_temp = condnum_all(indxWin+1, :);
    data_temp = log10(data_temp);
    [xCdf, yCdf] = func_cdf(data_temp(:));
    xCdf = 10.^xCdf;
    if indxWin < length(rollOffWin_set)
        semilogx(xCdf, yCdf, ['r-',marker_set{indxWin+1}], 'lineWidth', 1);
        hold on;
        rollOffWin = rollOffWin_set(indxWin+1);
        label_str = {['(proposed) Rx. RC Win.\&Overlap, $\alpha_{\mathrm{W}}=$', num2str(rollOffWin)]};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)
        semilogx(xCdf, yCdf, 'b--+', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 70dB-sidelobes'};
        label_set = [label_set, label_str];
    elseif indxWin == length(rollOffWin_set)+1
        semilogx(xCdf, yCdf, 'b--o', 'lineWidth', 1);
        hold on;
        label_str = {'(legacy) Tx. Chebyshev Win., 90dB-sidelobes'};
        label_set = [label_set, label_str];
    end
end
xlabel('Condition Number of $\mathbf{H}$', 'fontsize', 12, 'interpreter', 'latex');
ylabel('CDF', 'fontsize', 12, 'interpreter', 'latex');
legend(label_set, 'fontsize', 10, 'interpreter', 'latex', 'Box', 'Off');
grid on;
ax = gca;
ax.TickLabelInterpreter = 'latex';
axis([-inf,inf,-0.02,1.02]);
set(gcf, 'PaperType', 'a3');





%% AFDM transceiver
function [yData, yNoise] = afdmTransceiver(xData, a_wgn, ...
    numSc, lenCp, lenCp_ext, winTx, winRx, ...
    chirpParam1, chirpParam2, ...
    rollOffTime, filterOrder, upSampleCoef, ...
    numPath, tapGainPath, tapDelayPath, tapDopplerPath)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    %
    % time domain shaping pulse
    filterCoef = rootRaisedCosFilter(rollOffTime, filterOrder, upSampleCoef);
    filterCoefTx = filterCoef * sqrt(upSampleCoef);
    filterCoefRx = filterCoef / sqrt(upSampleCoef);
    % chirp base function
    chirpBase1 = exp(1i * 2 * pi * chirpParam1 * [-(lenCp+lenCp_ext): numSc-1].^2);
    chirpBase2 = exp(1i * 2 * pi * chirpParam2 * [0: numSc-1].^2);
    % prechirping and ifft
    s_ofdm = ifft(xData .* chirpBase2) * sqrt(numSc);
    % transmit shaping windowing
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    s_ofdmW = s_ofdm .* winTx;
    % cyclic-prefix attachment
    s_ofdmCp = [s_ofdmW(end-lenCp-lenCp_ext+1:end), s_ofdmW];
    % chirping
    s_afdm = s_ofdmCp .* chirpBase1;
    % transmission time domain pulse shaping
    numSample = numSc + lenCp + lenCp_ext;
    numSampleUp = upSampleCoef * numSample;
    s_t_upSamp = zeros(1, numSampleUp);
    s_t_upSamp(1: upSampleCoef: numSampleUp) = s_afdm;
    s_t_upSamp = conv(s_t_upSamp, filterCoefTx);
    % signals go through the doubly-selective channel
    s_r_upSamp = zeros(1, numSampleUp+filterOrder);
    for indxPath = 0: numPath-1
        % delay modulation
        phaseDelay = exp(-1i*2*pi*tapDelayPath(indxPath+1) ...
            * ifftshift([-(numSampleUp+filterOrder)/2:(numSampleUp+filterOrder)/2-1]) ...
            *upSampleCoef/(numSampleUp+filterOrder));
        s_pass = ifft(fft(s_t_upSamp) .* phaseDelay);
        % doppler modulation
        phaseDoppler = exp(1i*2*pi*tapDopplerPath(indxPath+1)/numSc ...
            * [-filterOrder/2-(lenCp+lenCp_ext)*upSampleCoef : ...
            numSampleUp-(lenCp+lenCp_ext)*upSampleCoef+filterOrder/2-1]/upSampleCoef);
        s_pass = s_pass .* phaseDoppler;
        % gain modulation
        s_r_upSamp = s_r_upSamp + tapGainPath(indxPath+1) * s_pass;
    end
    % receive time-domain pulse shaping
    s_r_upSamp = conv(s_r_upSamp, filterCoefRx); 
    s_r = s_r_upSamp(filterOrder+3: upSampleCoef: filterOrder+2+numSampleUp);
    % de-chirping
    s_r = s_r .* conj(chirpBase1);
    n_r = a_wgn .* conj(chirpBase1);
    % prefix removal
    s_r = s_r(lenCp+1: lenCp+lenCp_ext+numSc);
    n_r = n_r(lenCp+1: lenCp+lenCp_ext+numSc);
    % shaping windowing
    winRx = winRx / max(winRx);
    s_r = s_r .* winRx;
    n_r = n_r .* winRx;
    % time-domain overlap-summation
    s_r = s_r(lenCp_ext+1: end) + [zeros(1, numSc-lenCp_ext), s_r(1: lenCp_ext)];
    n_r = n_r(lenCp_ext+1: end) + [zeros(1, numSc-lenCp_ext), n_r(1: lenCp_ext)];
    % fft
    yData  = fft(s_r) / sqrt(numSc);
    yNoise = fft(n_r) / sqrt(numSc);
    % de-prechirping
    yData  = yData  .* conj(chirpBase2);
    yNoise = yNoise .* conj(chirpBase2);
end

%% Channel matrix in modulation-symbol domain
function hData = inOutAfdm_rxWin(numSc, indxSc_schd, lenCp_ext, ...
    chirpParam1, chirpParam2, rollOffTime, lenPulseTime, winRx, ...
    numPath, tapGainPath, tapDelayPath, tapDopplerPath)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % winRx = winRx / max(winRx);
    % numSc_schd = length(indxSc_schd);
    % hData = zeros(numSc, numSc_schd);
    % for indxPath = 0: numPath-1
    %     tapGain = tapGainPath(indxPath+1);
    %     tapDelay = tapDelayPath(indxPath+1);
    %     tapDoppler = tapDopplerPath(indxPath+1);
    %     for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
    %         pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
    %         for indxSc = 0: numSc_schd-1
    %             m_prime = indxSc_schd(indxSc+1);
    %             m_shift = m_prime - (2*numSc*chirpParam1*ell_prime - tapDoppler);
    %             winRx_phase = winRx ...
    %                 .* exp(1i * 2*pi/numSc * m_shift*[-lenCp_ext:numSc-1]);
    %             winRx_overlap = winRx_phase(lenCp_ext+1: end) ...
    %                 + [zeros(1, numSc-lenCp_ext), winRx_phase(1: lenCp_ext)];
    %             pulseCoefWin = fft(winRx_overlap) / numSc;
    %             chanCoef = tapGain ...
    %                 * exp(1i * 2*pi * chirpParam1*ell_prime^2) ...
    %                 * exp(- 1i * 2*pi/numSc * m_prime*ell_prime)...
    %                 * exp(1i * 2*pi * chirpParam2*(m_prime^2-[0:numSc-1].^2));
    %             hData(:, indxSc+1) = hData(:, indxSc+1) ...
    %                 + (chanCoef * pulseCoefTime .* pulseCoefWin).';
    %         end
    %     end
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    winRx = winRx / max(winRx);
    numSc_schd = length(indxSc_schd);
    hData = zeros(numSc+lenCp_ext, numSc_schd);
    for indxPath = 0: numPath-1
        tapGain = tapGainPath(indxPath+1);
        tapDelay = tapDelayPath(indxPath+1);
        tapDoppler = tapDopplerPath(indxPath+1);
        for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
            pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
            for indxSc = 0: numSc_schd-1
                m_prime = indxSc_schd(indxSc+1);
                m_shift = m_prime - (2*numSc*chirpParam1*ell_prime - tapDoppler);
                shiftPhase = exp(1i * 2*pi/numSc * m_shift*[-lenCp_ext:numSc-1]);
                chanCoef = tapGain ...
                    * exp(1i * 2*pi * chirpParam1*ell_prime^2) ...
                    * exp(- 1i * 2*pi/numSc * m_prime*ell_prime)...
                    * exp(1i * 2*pi * chirpParam2*m_prime^2) ...
                    / sqrt(numSc);
                hData(:, indxSc+1) = hData(:, indxSc+1) ...
                    + (chanCoef * pulseCoefTime .* shiftPhase).';
            end
        end
    end
    winRx = winRx.';
    for indxSc = 0: numSc_schd-1
        hData(:, indxSc+1) = winRx .* hData(:, indxSc+1);
    end
    hData = [hData(lenCp_ext+1: end, :)] + ...
        [zeros(numSc-lenCp_ext, numSc_schd); hData(1:lenCp_ext, :)];
    dePrechirpPhase = exp(- 1i * 2*pi * chirpParam2*[0:numSc-1].^2).';
    for indxSc = 0: numSc_schd-1
        hData(:, indxSc+1) = fft(hData(:, indxSc+1)) / sqrt(numSc);
        hData(:, indxSc+1) = dePrechirpPhase .* hData(:, indxSc+1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % numSc_schd = length(indxSc_schd);
    % lenPulseWin = 32;
    % hData = zeros(numSc, numSc_schd);
    % for indxPath = 0: numPath-1
    %     tapGain = tapGainPath(indxPath+1);
    %     tapDelay = tapDelayPath(indxPath+1);
    %     tapDoppler = tapDopplerPath(indxPath+1);
    %     for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
    %         pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
    %         for indxSc = 0: numSc_schd-1
    %             m_prime = indxSc_schd(indxSc+1);
    %             m_shift = m_prime - (2*numSc*chirpParam1*ell_prime - tapDoppler);
    %             for m = floor(m_shift)-lenPulseWin/2: ceil(m_shift)+lenPulseWin/2
    %                 pulseCoefWin = raisedCosPulse(rollOffWin, m-m_shift);
    %                 chanCoef = tapGain ...
    %                     * exp(- 1i * 2*pi/numSc * ((numSc-lenCp_ext)/2)*(m-m_shift)) ...
    %                     * exp(1i * 2*pi * chirpParam1*ell_prime^2) ...
    %                     * exp(- 1i * 2*pi/numSc * m_prime*ell_prime)...
    %                     * exp(1i * 2*pi * chirpParam2*(m_prime^2-m.^2));
    %                 hData(mod(m,numSc)+1, indxSc+1) = hData(mod(m,numSc)+1, indxSc+1) ...
    %                     + (chanCoef * pulseCoefTime .* pulseCoefWin).';
    %             end
    %         end
    %     end
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function hData = inOutAfdm_txWin(numSc, indxSc_schd, ...
    chirpParam1, chirpParam2, rollOffTime, lenPulseTime, winTx, ...
    numPath, tapGainPath, tapDelayPath, tapDopplerPath)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % winTx = winTx / sqrt(mean(abs(winTx).^2));
    % numSc_schd = length(indxSc_schd);
    % hData = zeros(numSc, numSc_schd);
    % for indxPath = 0: numPath-1
    %     tapGain = tapGainPath(indxPath+1);
    %     tapDelay = tapDelayPath(indxPath+1);
    %     tapDoppler = tapDopplerPath(indxPath+1);
    %     for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
    %         pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
    %         for indxSc = 0: numSc_schd-1
    %             m_prime = indxSc_schd(indxSc+1);
    %             m_shift = m_prime - (2*numSc*chirpParam1*ell_prime - tapDoppler);
    %             winTx_shift = circshift(winTx, ell_prime);
    %             winTx_phase = winTx_shift .* exp(1i * 2*pi/numSc * m_shift*[0:numSc-1]);
    %             pulseCoefWin = fft(winTx_phase) / numSc;
    %             chanCoef = tapGain ...
    %                 * exp(1i * 2*pi * chirpParam1*ell_prime^2) ...
    %                 * exp(- 1i * 2*pi/numSc * m_prime*ell_prime)...
    %                 * exp(1i * 2*pi * chirpParam2*(m_prime^2-[0:numSc-1].^2));
    %             hData(:, indxSc+1) = hData(:, indxSc+1) ...
    %                 + (chanCoef * pulseCoefTime .* pulseCoefWin).';
    %         end
    %     end
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    numSc_schd = length(indxSc_schd);
    hData = zeros(numSc, numSc_schd);
    for indxPath = 0: numPath-1
        tapGain = tapGainPath(indxPath+1);
        tapDelay = tapDelayPath(indxPath+1);
        tapDoppler = tapDopplerPath(indxPath+1);
        for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
            pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
            for indxSc = 0: numSc_schd-1
                m_prime = indxSc_schd(indxSc+1);
                m_shift = m_prime - (2*numSc*chirpParam1*ell_prime - tapDoppler);
                winTx_shift = circshift(winTx, ell_prime);
                winTx_phase = winTx_shift .* exp(1i * 2*pi/numSc * m_shift*[0:numSc-1]);
                chanCoef = tapGain ...
                    * exp(1i * 2*pi * chirpParam1*ell_prime^2) ...
                    * exp(- 1i * 2*pi/numSc * m_prime*ell_prime)...
                    * exp(1i * 2*pi * chirpParam2*m_prime^2) ...
                    / sqrt(numSc);
                hData(:, indxSc+1) = hData(:, indxSc+1) ...
                    + (chanCoef * pulseCoefTime .* winTx_phase).';
            end
        end
    end
    dePrechirpPhase = exp(- 1i * 2*pi * chirpParam2*[0:numSc-1].^2).';
    for indxSc = 0: numSc_schd-1
        hData(:, indxSc+1) = fft(hData(:, indxSc+1)) / sqrt(numSc);
        hData(:, indxSc+1) = dePrechirpPhase .* hData(:, indxSc+1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Channel estimation using embedded pilot
function hMat_est = chanMat_est(yPilot, xPilot, indxPilot, ...
    numSc, indxSc_schd, chirpParam1, chirpParam2, labelDelay, labelDoppler)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    numDoppler = length(labelDoppler);
    numDelay = length(labelDelay);
    hDD_est = zeros(numDoppler, numDelay);
    for indxDoppler = 0: numDoppler-1
        for indxDelay = 0: numDelay-1
            valDoppler = labelDoppler(indxDoppler+1);
            valDelay = labelDelay(indxDelay+1);
            indxX_shift = 2*numSc*chirpParam1*valDelay - valDoppler;
            indxY = mod(indxPilot-indxX_shift, numSc);
            phaseCompTerm = exp(1i * 2*pi * chirpParam2 * (indxPilot^2 - indxY.^2)) ...
                * exp(1i * 2*pi * chirpParam1*valDelay^2) ...
                * exp(- 1i * 2*pi/numSc * indxPilot*valDelay);
            hDD_est(indxDoppler+1, indxDelay+1) ...
                = yPilot(indxY+1) / xPilot(indxPilot+1) * conj(phaseCompTerm);
        end
    end
    numSc_schd = length(indxSc_schd);
    hMat_est = zeros(numSc, numSc_schd);
    for indxDelay = 0: numDelay-1
        for indxDoppler = 0: numDoppler-1
            tapGain = hDD_est(indxDoppler+1, indxDelay+1);
            tapDelay = labelDelay(indxDelay+1);
            tapDoppler = labelDoppler(indxDoppler+1);
            indxX_shift = 2*numSc*chirpParam1*tapDelay - tapDoppler;
            chanCoef = tapGain ...
                    * exp(1i * 2*pi/numSc * tapDelay*tapDoppler) ...
                    * exp(- 1i * 2*pi * chirpParam1*tapDelay^2) ...
                    * exp(- 1i * 2*pi * chirpParam2*indxX_shift^2);
            for indxSc = 0: numSc_schd-1
                indxX = indxSc_schd(indxSc+1);
                indxY = mod(indxX - indxX_shift, numSc);
                phaseTerm = exp(1i * 2*pi/numSc * (2*numSc*chirpParam2*indxX_shift-tapDelay)*indxY);
                hMat_est(indxY+1, indxSc+1) = hMat_est(indxY+1, indxSc+1) ...
                    + chanCoef * phaseTerm;
            end
        end
    end
    % figure; bar3(abs(hDD_est).^2);
    % xlabel('Delay bins'); ylabel('Doppler bins'); zlabel('Magnitude');
    % set(gca, 'xtick',[1:4:numDelay], 'xticklabel',labelDelay(1:4:end));
    % set(gca, 'ytick',[1:numDoppler], 'yticklabel',labelDoppler);
end
    
%% Root raised cosine pulse shaping filter
function filterCoef = rootRaisedCosFilter(rollOff, filterOrder, upSampleCoef)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    indxFreq = [-filterOrder/2: filterOrder/2] / (filterOrder+1);
    winFilter = zeros(1, filterOrder+1);
    winFilter(abs(indxFreq) <= (1-rollOff)/(2*upSampleCoef)) = 1;
    indxFreq_other = indxFreq((abs(indxFreq)>(1-rollOff)/(2*upSampleCoef)) ...
        &(abs(indxFreq)<=(1+rollOff)/(2*upSampleCoef)));
    winFilter((abs(indxFreq)>(1-rollOff)/(2*upSampleCoef)) ...
        &(abs(indxFreq)<=(1+rollOff)/(2*upSampleCoef))) ...
        = cos(pi/(2*rollOff)*(abs(indxFreq_other)*upSampleCoef-(1-rollOff)/2)).^2;
    winFilter = sqrt(winFilter);
    filterCoef = ifftshift(ifft(ifftshift(winFilter)));
    filterCoef = filterCoef / sqrt(sum(abs(filterCoef).^2));
end

%% Raised cosine pulse shaping window
function winRx = raisedCosWin(rollOff, numSc)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    indxTime = [-(1+rollOff)*numSc/2: (1+rollOff)*numSc/2-1];
    winRx = zeros(1, (1+rollOff)*numSc);
    winRx(abs(indxTime) <= (1-rollOff)*numSc/2) = 1;
    indxTime_other = indxTime(~(abs(indxTime) <= (1-rollOff)*numSc/2));
    winRx(~ (abs(indxTime) <= (1-rollOff)*numSc/2)) ...
        = cos(pi/(2*rollOff*numSc)*(abs(indxTime_other)-(1-rollOff)*numSc/2)).^2;
    winRx = winRx / sqrt(mean(abs(winRx).^2));
end

%% Raised cosine pulse
function pulseCoef = raisedCosPulse(rollOff, indx)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    term_0 = sin(pi*indx) ./ (pi*indx);
    term_0(indx==0) = 1;
    term_1 = cos(rollOff*pi*indx) ./ (1 - 4*rollOff^2*indx.^2 + eps);
	term_1(abs(2*rollOff*indx)==1) = pi / 4;
    pulseCoef = term_0 .* term_1;
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
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    modulatOrder = length(bits);  % scheduled modulation order
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

%% CDF
function [xCdf, yCdf] = func_cdf(data)
    % Coded by Haojian Zhang
    % UWB-LAB, School of Information Science and Technology, Harbin Institute of Technology, Shenzhen
    % Copyright (c) 2025, all rights reserved.
    numData = numel(data);
    numCdf = ceil(sqrt(numData));
    xCdf = linspace(min(data), max(data), numCdf);
    yCdf = zeros(1, numCdf);
    for indxX = 0: numCdf-1
        yCdf(indxX+1) = sum(data <= xCdf(indxX+1)) / numData;
    end
end