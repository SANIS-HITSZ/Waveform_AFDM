% Simulation of the legacy PS-AFDM transceiver
% Coded by Haojian Zhang
% UWB-LAB, Rehool of Information Reience and Technology, Harbin Institute of Technology, Shenzhen
% Copyright (c) 2025, all rights reserved.
clc; clear; close all;
% rng(1013);

%% Parameters of signalling
numRe = 4096;         % number of resource elements in resource grid
lenCp = 288;          % length of normal cyclic-prefix
modulatOrder = 10;    % scheduled modulation order
filterOrder = 128;   % order of the time-domain pulse shaping filter (after up-sampling)
rollOffTime = 0.2;    % roll-off factor for the time-domain pulse shaping (root-raised-cosine)
upSampleCoef = 4;     % upsampling factor in the time domain
lenPulseTime = 32;    % considered valid length of the overall time-domain shaping pulse
bitsRef_temp = zeros(2^modulatOrder, modulatOrder);    % transmit bit-alphabet
qamsRef_temp = zeros(1, 2^modulatOrder);               % transmit qam-alphabet
for indxQam = 0: 2^modulatOrder-1
    bitsRef_temp(indxQam+1, :) = de2bi(indxQam, modulatOrder, 'left-msb');
    qamsRef_temp(indxQam+1) = func_nrQamMapper(bitsRef_temp(indxQam+1, :));
end
startRe_schd = 0;     % scheduled start resource element index
numRe_schd = 50*12;   % scheduled resource element number
indxRe_schd = [startRe_schd: startRe_schd + numRe_schd - 1];  % scheduled resource element indexes
typeChanMat = 1;

%% Monte Carlo tests
numTest = 10000;                   % the number of monte carlo tests
snr_set = [20: 5: 55];             % received snr (db)
numPath_set = [3, 10];
ell_max_set = [1, 10, 100];
k_max_set = [1, 3, 6];
winTx_set = [0, 1, 2]; % ["hamm","cheb-70","cheb-90"];

params_set = [];
for i_sinr = 0: length(snr_set)-1
    for i_num = 0: length(numPath_set)-1
        for i_ell = 0: length(ell_max_set)-1
            for i_k = 0: length(k_max_set)-1
                for i_win = 0: length(winTx_set)-1
                    params_set = [ params_set, [snr_set(i_sinr+1); ...
                                                numPath_set(i_num+1); ...
                                                ell_max_set(i_ell+1); ...
                                                k_max_set(i_k+1); ...
                                                winTx_set(i_win+1)] ];
                end
            end
        end
    end
end

numParams = size(params_set, 2);


ber_all_tx = zeros(1, numParams);
nmse_all_tx = zeros(1, numParams);
numTest_all_tx = zeros(1, numParams);


% delete(gcp('nocreate'));
% ppooll = parpool(10);


% parfor indxParams = 0: numParams-1
for indxParams = 0: numParams-1

    params_vec = params_set(:, indxParams+1);
    snrDb = params_vec(1);
    numPath = params_vec(2);
    ell_max = params_vec(3);
    k_max = params_vec(4);
    winRx_type = params_vec(5);
    bitsRef = bitsRef_temp;
    qamsRef = qamsRef_temp;

    nmseEst = 0;
    berEst = 0;
    cntErr = 0;

    rng(1013);
    for indxTest = 0: numTest-1

        %% Parameters of delay-Doppler channel
        tapDelayPath = rand(1,numPath)*ell_max + lenPulseTime/2;  % delay taps
        tapDopplerPath = cos(rand(1,numPath)*2*pi)*k_max;         % doppler taps
        tapGainPath = (randn(1,numPath)+1i*randn(1,numPath))/sqrt(numPath*2);  % rayleigh fading path gains

        %% Chirp-rate
        k_reserve = 4;     % reserved guard spacing for fractional doppler
        chirpRate = (2*(ceil(k_max)+k_reserve)+1)/(2*numRe);  % chirp-rate
        prechirpRate = 0;   % prechirp-rate

        %% Data modulation symbols
        symData_user = round(rand(1, numRe_schd)*2^modulatOrder - 0.5);   % transmit decimal symbols
        bitData_user = bitsRef(symData_user+1, :);  % transmit bits
        xData_user   = qamsRef(symData_user+1);     % transmit qams

        %% Resource grid mapping
        numDoppler = abs(2*numRe*chirpRate);
        numDelay = ceil(ell_max) + lenPulseTime + 1;
        xAll = zeros(1, numRe);                 % resource grid
        % data
        indxRegionData = mod([min(indxRe_schd)-numDelay*numDoppler+(numDoppler-1)/2+1: max(indxRe_schd)+(numDoppler-1)/2], numRe);
        xAll(indxRe_schd+1) = xData_user;       % data mapping
        % pilot
        indxPilot = mod(min(indxRe_schd)-numDelay*numDoppler, numRe);
        indxRegionPilot = mod([indxPilot-numDelay*numDoppler+(numDoppler-1)/2+1: indxPilot+(numDoppler-1)/2], numRe);
        xPilot = sqrt(length(indxRegionData)+length(indxRegionPilot)-numRe_schd);
        xAll(indxPilot+1) = xPilot;                 % pilot mapping
        %
        if sum(ismember(indxRegionPilot,indxRegionData), 'all') > 0
            disp(['The channel dispersion is too large. ', ...
                num2str(sum(ismember(indxRegionPilot,indxRegionData), 'all')), ...
                '. The pilot and data will overlap after passing through the channel.']);
            disp('At this point, it is recommended to place the pilot in an AFDM symbol that is different from the data.');
            break;
        end

        %% Shaping window
        lenCp_ext = 0;   % length of extended cyclic-prefix for receiving windowing
        if winRx_type == 0
            winTx = hamming(numRe).';
        elseif winRx_type == 1
            winTx = chebwin(numRe, 70).';   % transmit shaping window
        else
            winTx = chebwin(numRe, 90).';
        end
        winRx = rectwin(numRe).';       % receive shaping window

        %% AFDM transceiver procedures
        numSample = numRe + lenCp + lenCp_ext;
        a_wgn = (randn(1,numSample)+1i*randn(1,numSample)) / sqrt(2);  % normalized additive white noise
        [yAll, yNoiseAll] = afdmTransceiver(xAll, a_wgn, ...
            numRe, lenCp, lenCp_ext, winTx, winRx, ...
            chirpRate, prechirpRate, ...
            rollOffTime, filterOrder, upSampleCoef, ...
            numPath, tapGainPath, tapDelayPath, tapDopplerPath);
        % pilot
        yPilot = yAll(indxRegionPilot+1);
        yNoisePilot = yNoiseAll(indxRegionPilot+1);
        yPilot_noised = yPilot + sqrt(10^(-snrDb/10))*yNoisePilot;
        powerNoisePilot = mean(abs(sqrt(10^(-snrDb/10))*yNoisePilot).^2);
        % data
        yData = yAll(indxRegionData+1);
        yNoiseData = yNoiseAll(indxRegionData+1);
        yData_noised = yData + sqrt(10^(-snrDb/10))*yNoiseData;
        powerNoiseData = mean(abs(sqrt(10^(-snrDb/10))*yNoiseData).^2);

        %% Channel matrix construction
        if typeChanMat == 0
            hData = chanMatAfdm(numRe, indxRe_schd, lenCp_ext, ...
                chirpRate, prechirpRate, rollOffTime, winTx, winRx, ...
                numPath, tapGainPath, tapDelayPath, tapDopplerPath);
        else
            numDoppler = abs(2 * numRe * chirpRate);
            numDelay = ceil(ell_max) + lenPulseTime + 1;
            labelDoppler = [-(numDoppler-1)/2: (numDoppler-1)/2];
            labelDelay = [0: numDelay-1];
            [hData, hDD_est] = chanMat_est(yPilot_noised, xPilot, ...
                indxPilot, indxRegionPilot, powerNoisePilot, ...
                numRe, indxRe_schd, chirpRate, prechirpRate, labelDelay, labelDoppler);
        end
        hData = hData(indxRegionData+1, :);

        %% LMMSE equalization
        xDataEst = (hData'*hData + powerNoiseData*eye(numRe_schd)) \ (hData'*yData_noised.');
        bitEst = zeros(numRe_schd, modulatOrder);
        for indxX = 0: numRe_schd-1
            [~, indxMin] = min(abs(xDataEst(indxX+1) - qamsRef));
            bitEst(indxX+1, :) = bitsRef(indxMin, :);
        end

        %%
        nmse_temp = sum(abs(yData.' - hData*xData_user.').^2) / sum(abs(yData).^2);
        errRate = sum(sum(bitEst~=bitData_user)) / (numRe_schd*modulatOrder);
        nmseEst = nmseEst + nmse_temp;
        berEst = berEst + errRate;
        cntErr = cntErr + double(errRate>0);
        if cntErr > 400
            break;
        end

    end

    nmse_all_tx(indxParams+1) = nmseEst / (indxTest+1);
    ber_all_tx(indxParams+1) = berEst / (indxTest+1);
    numTest_all_tx(indxParams+1) = indxTest+1;

    disp([num2str(indxParams+1),' / ',num2str(numParams)]);

end

save('data_ber_all_tx.mat', 'ber_all_tx');
save('data_nmse_all_tx.mat', 'nmse_all_tx');
save('data_numTest_all_tx.mat', 'numTest_all_tx');

delete(gcp('nocreate'));





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% Some Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AFDM transceiver
function [yAll, yNoiseAll] = afdmTransceiver(xAll, a_wgn, ...
        numRe, lenCp, lenCp_ext, winTx, winRx, ...
        chirpRate, prechirpRate, ...
        rollOffTime, filterOrder, upSampleCoef, ...
        numPath, tapGainPath, tapDelayPath, tapDopplerPath)
    % time domain shaping pulse
    filterCoef = rootRaisedCosFilter(rollOffTime, filterOrder, upSampleCoef);
    filterCoefTx = filterCoef * sqrt(upSampleCoef);
    filterCoefRx = filterCoef / sqrt(upSampleCoef);
    % chirp base function
    chirpBase1 = exp(1i * 2 * pi * chirpRate * [-(lenCp+lenCp_ext): numRe-1].^2);
    chirpBase2 = exp(1i * 2 * pi * prechirpRate * [0: numRe-1].^2);
    % prechirping and ifft
    s_ofdm = ifft(xAll .* chirpBase2) * sqrt(numRe);
    % transmit shaping windowing
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    s_ofdmW = s_ofdm .* winTx;
    % cyclic-prefix attachment
    s_ofdmCp = [s_ofdmW(end-lenCp-lenCp_ext+1:end), s_ofdmW];
    % chirping
    s_afdm = s_ofdmCp .* chirpBase1;
    % transmission time domain pulse shaping
    numSample = numRe + lenCp + lenCp_ext;
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
        phaseDoppler = exp(1i*2*pi*tapDopplerPath(indxPath+1)/numRe ...
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
    s_r = s_r(lenCp+1: lenCp+lenCp_ext+numRe);
    n_r = n_r(lenCp+1: lenCp+lenCp_ext+numRe);
    % receive shaping windowing
    winRx = winRx / max(abs(winRx));
    s_r = s_r .* winRx;
    n_r = n_r .* winRx;
    % time-domain overlap-summation
    s_r = s_r(lenCp_ext+1: end) + [zeros(1, numRe-lenCp_ext), s_r(1: lenCp_ext)];
    n_r = n_r(lenCp_ext+1: end) + [zeros(1, numRe-lenCp_ext), n_r(1: lenCp_ext)];
    % fft
    yAll  = fft(s_r) / sqrt(numRe);
    yNoiseAll = fft(n_r) / sqrt(numRe);
    % de-prechirping
    yAll      = yAll      .* conj(chirpBase2);
    yNoiseAll = yNoiseAll .* conj(chirpBase2);
end

%% Channel matrix in modulation-symbol domain
function hData = chanMatAfdm(numRe, indxRe_schd, lenCp_ext, ...
        chirpRate, prechirpRate, rollOffTime, winTx, winRx, ...
        numPath, tapGainPath, tapDelayPath, tapDopplerPath)
    lenPulseTime = 32;
    winTx = winTx / sqrt(mean(abs(winTx).^2));
    winTx = [winTx(end-lenCp_ext+1:end), winTx];
    winRx = winRx.' / max(winRx);
    numRe_schd = length(indxRe_schd);
    hData = zeros(numRe+lenCp_ext, numRe_schd);
    for indxPath = 0: numPath-1
        tapGain = tapGainPath(indxPath+1);
        tapDelay = tapDelayPath(indxPath+1);
        tapDoppler = tapDopplerPath(indxPath+1);
        for ell_prime = floor(tapDelay)-lenPulseTime/2: ceil(tapDelay)+lenPulseTime/2
            pulseCoefTime = raisedCosPulse(rollOffTime, ell_prime-tapDelay);
            for indxRe = 0: numRe_schd-1
                m_prime = indxRe_schd(indxRe+1);
                m_shift = m_prime - (2*numRe*chirpRate*ell_prime - tapDoppler);
                winTx_shift = circshift(winTx, ell_prime);
                winTx_phase = winTx_shift .* exp(1i * 2*pi/numRe * m_shift*[-lenCp_ext:numRe-1]);
                chanCoef = tapGain ...
                    * exp(1i * 2*pi * chirpRate*ell_prime^2) ...
                    * exp(- 1i * 2*pi/numRe * m_prime*ell_prime)...
                    * exp(1i * 2*pi * prechirpRate*m_prime^2) ...
                    / sqrt(numRe);
                hData(:, indxRe+1) = hData(:, indxRe+1) ...
                    + (chanCoef * pulseCoefTime .* winTx_phase).';
            end
        end
    end
    hData = winRx .* hData;
    hData = [hData(lenCp_ext+1: end, :)] + ...
        [zeros(numRe-lenCp_ext, numRe_schd); hData(1:lenCp_ext, :)];
    dePrechirpPhase = exp(- 1i * 2*pi * prechirpRate*[0:numRe-1].^2).';
	hData = fft(hData, [], 1) / sqrt(numRe);
	hData = dePrechirpPhase .* hData;
end

%% Channel estimation using embedded pilot
function [hMat_est, hDD_est] = chanMat_est(yPilot, xPilot, ...
    indxPilot, indxRegion, powerNoise, ...
    numRe, indxRe_schd, chirpRate, prechirpRate, labelDelay, labelDoppler)
    % threshold filtering
    probFalseAlarm = 1e-4;
    thres = sqrt( - powerNoise*log(probFalseAlarm) );
    yPilot_filt = yPilot;
    yPilot_filt(abs(yPilot) <= thres) = 0;
    yPilot = zeros(1, numRe);
    yPilot(indxRegion+1) = yPilot_filt;
    % delay-doppler channel
    numDoppler = length(labelDoppler);
    numDelay = length(labelDelay);
    hDD_est = zeros(numDoppler, numDelay);
    for indxDelay = 0: numDelay-1
        for indxDoppler = 0: numDoppler-1
            valDoppler = labelDoppler(indxDoppler+1);
            valDelay = labelDelay(indxDelay+1);
            indxX_shift = 2*numRe*chirpRate*valDelay - valDoppler;
            indxY = mod(indxPilot-indxX_shift, numRe);
            phaseCompTerm = exp(1i * 2*pi * chirpRate*valDelay^2) ...
                * exp(1i * 2*pi * prechirpRate * (indxPilot^2-indxY.^2)) ...
                * exp(- 1i * 2*pi/numRe * indxPilot*valDelay);
            hDD_est(indxDoppler+1, indxDelay+1) ...
                = yPilot(indxY+1) / xPilot * conj(phaseCompTerm);
        end
    end
    % channel matrix
    numRe_schd = length(indxRe_schd);
    hMat_est = zeros(numRe, numRe_schd);
    for indxDelay = 0: numDelay-1
        for indxDoppler = 0: numDoppler-1
            tapGain = hDD_est(indxDoppler+1, indxDelay+1);
            tapDelay = labelDelay(indxDelay+1);
            tapDoppler = labelDoppler(indxDoppler+1);
            indxX_shift = 2*numRe*chirpRate*tapDelay - tapDoppler;
            for indxRe = 0: numRe_schd-1
                indxX = indxRe_schd(indxRe+1);
                indxY = mod(indxX - indxX_shift, numRe);
                phaseTerm = exp(1i * 2*pi * chirpRate*tapDelay^2) ...
                    * exp(1i * 2*pi * prechirpRate * (indxX^2-indxY^2)) ...
                    * exp(-1i * 2*pi/numRe * indxX*tapDelay);
                hMat_est(indxY+1, indxRe+1) = hMat_est(indxY+1, indxRe+1) ...
                    + tapGain * phaseTerm;
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
function winRx = raisedCosWin(rollOff, numRe)
    indxTime = [-(1+rollOff)*numRe/2: (1+rollOff)*numRe/2-1];
    winRx = zeros(1, (1+rollOff)*numRe);
    winRx(abs(indxTime) <= (1-rollOff)*numRe/2) = 1;
    indxTime_other = indxTime(~(abs(indxTime) <= (1-rollOff)*numRe/2));
    winRx(~ (abs(indxTime) <= (1-rollOff)*numRe/2)) ...
        = cos(pi/(2*rollOff*numRe)*(abs(indxTime_other)-(1-rollOff)*numRe/2)).^2;
    winRx = winRx / sqrt(mean(abs(winRx).^2));
end

%% Raised cosine pulse
function pulseCoef = raisedCosPulse(rollOff, indx)
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
    % UWB-LAB, Rehool of Information Reience and Technology, Harbin Institute of Technology, Shenzhen
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

%%  Channel-awared receive pulse shaping window
function [numTap, powerTap, delayTap, dopplerTap] ...
    = chanTap_aware(yPilot, powerNoise, indxRegionPilot, ...
        prechirpRate, upSampCoefDaft, numDelay, numDoppler)
    % threshold filtering
    probFalseAlarm = 1e-4;
    thres = sqrt( - powerNoise*log(probFalseAlarm) );
    yPilot_filt = yPilot;
    yPilot_filt(abs(yPilot) <= thres) = 0;
    % upsampling
    yPilot_upSamp = fft(ifft(yPilot_filt).*chebwin(numDelay*numDoppler,70).', ...
        upSampCoefDaft*numDelay*numDoppler);
    % peaks extraction
    [pks, locs] = findpeaks(abs(yPilot_upSamp));
    indxDoppler = mod(locs-1, upSampCoefDaft*numDoppler);
    indxDelay = floor((locs-1)/(upSampCoefDaft*numDoppler));
    labelDoppler = [-upSampCoefDaft*(numDoppler-1)/2:upSampCoefDaft*(numDoppler+1)/2-1]/upSampCoefDaft;
    labelDelay = fliplr([0: numDelay-1]);
    powerTap = abs(pks).^2;
    delayTap = labelDelay(indxDelay+1);
    dopplerTap = labelDoppler(indxDoppler+1);
    numTap = length(powerTap);
end


function winRx = chanAwareWolaWinRx(winRx, numRe, lenCp_ext, chirpRate, ...
    winOrder, numPath, powerPath, dopplerPath)
    % base functions
    winRx = winRx.';
    cosBase = zeros(lenCp_ext, winOrder);
    indxBase = [0:lenCp_ext-1].' / (lenCp_ext-1);
    for indxOrder = 0: winOrder-1
        cosBase(:,indxOrder+1) = indxBase .* (1-indxBase) ...
            .* cos(indxOrder*pi*indxBase);
    end
    winBase_prefix_0 = winRx(1:lenCp_ext);
    winBase_prefix_1 = cosBase;
    winBase_main = winRx(lenCp_ext+1:end-lenCp_ext);
    winBase_tail_0 = winRx(end-lenCp_ext+1:end);
    winBase_tail_1 = cosBase;
    % fourier operators
    numDoppler = abs(2 * numRe * chirpRate);
    operatorParFourierInn_main = exp(-1i*2*pi/numRe ...
        * [-(numDoppler-1)/2:(numDoppler-1)/2].'*[0:numRe-lenCp_ext-1]) / sqrt(numRe);
    operatorParFourierInn_tail = exp(-1i*2*pi/numRe ...
        * [-(numDoppler-1)/2:(numDoppler-1)/2].'*[numRe-lenCp_ext:numRe-1]) / sqrt(numRe);
    % energy matrix
    rInn = 0;
    gram_00 = 0;
    gram_11 = 0;
    gram_22 = 0;
    gram_01 = 0;
    gram_02 = 0;
    gram_12 = 0;
    for indxPath = 0: numPath-1
        valPower = powerPath(indxPath+1);
        valDoppler = dopplerPath(indxPath+1);
        dopplerCoef_prefix = exp(1i*2*pi/numRe*valDoppler*[-lenCp_ext:-1].');
        dopplerCoef_main = exp(1i*2*pi/numRe*valDoppler*[0:numRe-lenCp_ext-1].');
        dopplerCoef_tail = exp(1i*2*pi/numRe*valDoppler*[numRe-lenCp_ext:numRe-1].');
        winDoppler_prefix_0 = dopplerCoef_prefix .* winBase_prefix_0;
        winDoppler_prefix_1 = dopplerCoef_prefix .* winBase_prefix_1;
        winDoppler_main = dopplerCoef_main .* winBase_main;
        winDoppler_tail_0 = dopplerCoef_tail .* winBase_tail_0;
        winDoppler_tail_1 = dopplerCoef_tail .* winBase_tail_1;
        winDoppler_ola_00 = winDoppler_main;
        winDoppler_ola_10 = winDoppler_prefix_0 + winDoppler_tail_0;
        winDoppler_ola_11 = winDoppler_prefix_1;
        winDoppler_ola_12 = winDoppler_tail_1;
        winDftInn = [operatorParFourierInn_main*winDoppler_ola_00 ...
            + operatorParFourierInn_tail*winDoppler_ola_10, ...
            operatorParFourierInn_tail*winDoppler_ola_11, ...
            operatorParFourierInn_tail*winDoppler_ola_12];
        rInn = rInn + valPower*(winDftInn'*winDftInn);
        gram_00 = gram_00 + valPower*(winDoppler_ola_00'*winDoppler_ola_00 ...
                                    + winDoppler_ola_10'*winDoppler_ola_10);
        gram_11 = gram_11 + valPower*(winDoppler_ola_11'*winDoppler_ola_11);
        gram_22 = gram_22 + valPower*(winDoppler_ola_12'*winDoppler_ola_12);
        gram_01 = gram_01 + valPower*(winDoppler_ola_10'*winDoppler_ola_11);
        gram_02 = gram_02 + valPower*(winDoppler_ola_10'*winDoppler_ola_12);
        gram_12 = gram_12 + valPower*(winDoppler_ola_11'*winDoppler_ola_12);
    end
    rAll = [gram_00,  gram_01,  gram_02; ...
            gram_01', gram_11,  gram_12; ...
            gram_02', gram_12', gram_22];
    rAll = rAll + 1e-5*mean(abs(diag(rAll)))*eye(2*winOrder+1);
    % window coefficient solutions
    [weightMatrix, weightValue] = eig(rInn, rAll);
    [~, indxWeight] = max(real(diag(weightValue)));
    weightCoef = weightMatrix(:, indxWeight(1));
    weightCoef = weightCoef / weightCoef(1);
    winRx = [winBase_prefix_0 + winBase_prefix_1*weightCoef(2:winOrder+1); ...
             winBase_main; ...
             winBase_tail_0 + winBase_tail_1*weightCoef(winOrder+2:end)].';
end
