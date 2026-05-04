function [SpectData] = calibrateONChydrophoneSpectra(Y,Fs,S,df,overlapFactor,specType,Tavg)
    %---------------------------------------------------------------------
    % calibrateONChydrophoneSpectra.m
    % Description: Process Ocean Networks Canada single-channel hydrophone audio data to get calibrated power spectrum or power spectral density
    %
    % Syntax:
    % [SpectData] = calibrateONChydrophoneSpectra(audioPath,calPath,overlapFactor,specType)
    %
    % Input:    Y: audio time series [counts] - note this must be scaled to digital counts for the result to be calibrated properly
    %           Fs: audio sample rate [Hz]
    %           S: sensitivity data, formatted as a matrix with 1st column: frequency [Hz], 2nd column: sensitivity value [dB re counts/uPa]
    %           df: desired frequency bin spacing [Hz]. Sensitivity file will be interpolated if necessary.
    %           overlapFactor: ratio of overlap between adjacent time windows, typically 0 or 0.5
    %           specType: either 'power' for power spectrum [dB re 1µPa^2] or 'psd' for power spectral densitiy [dB re 1µPa^2/Hz]
    %           Tavg: time averaging duration [s]
    %
    % Output:   SpectData structure containing
    %           T: Time vector [s] corresponding to time steps in P
    %           F: Frequency vector [Hz]
    %           P: calibrated power spectrum [dB re 1µPa^2] or power spectral density [dB re 1µPa^2/Hz] depending on specified specType
    %           Pavg: time-average over all windows of P, units same as P
    %           PavgT: time-average over duration specified by Tavg, units same as P
    %           specType: 'power spectrum' or 'power spectral density' depending on specified specType
    %           specUnits: units of P and Pavg, either 'dB re 1µPa^2' for power spectrum or 'dB re 1µPa^2/Hz' for power spectral density, depending on specified specType
    %
    % References:   Merchant et al. (2015) Measuring acoustic habitats, Appendix S1. (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12330)
    %               Cerna, M. & Harvey, A.F. (2000) The fundamentals of FFT-based signal analysis and measurements, Application Note 041. Tech. rep. (https://www.sjsu.edu/people/burford.furman/docs/me120/FFT_tutorial_NI.pdf)
    %
    % Date/Author:  2026-04-23, Brendan Smith
    %---------------------------------------------------------------------

    % interpolate sensitivity curve, if df is different than existing frequency resolution
    Slin = 10.^(S(:,2)./20); % convert sensitivity values to linear space for interpolation [counts/µPa]
    Ftarget = S(1,1):df:S(end,1); % new sensitivity frequency vector based on desired spectral resolution [Hz]
    Sinterp_lin = interp1(S(:,1),Slin,Ftarget); % interpolated sensitivity values in linear space [counts/µPa]
    Sinterp(:,1) = Ftarget; % assign interpolated frequency vector to new sensitivity array [Hz]
    Sinterp(:,2) = 20*log10(Sinterp_lin); % assign interpolated sensitivity values (converted to dB) to new sensitivity array [dB re counts/µPa]
    S = Sinterp; % overwrite S with interpolated version

    % calculate spectral analysis parameters and define constants
    NFFT = round(Fs/df); % FFT size [samples]
    windowFn = hann(NFFT,'periodic'); % window function
    alpha = 0.5; % coherent gain factor (scaling factor) of hann window, which corrects for the reduction in amplitude introduced by the window function (Cerna & Harvey, 2000)
    
    % divide signal into time segments
    Ywin = buffer(Y,NFFT,round(overlapFactor*NFFT));
    Ywin = Ywin .* windowFn / alpha; % multiply each segment by the window function, and divide by window scaling factor
    
    % compute power spectrum
    X = fft(Ywin,NFFT,1); % discrete Fourier transform (DFT) of windowed data, operating on each column
    P = abs(X/NFFT).^2; % power spectrum
    Pss = P(1:floor(NFFT/2),:); % single-sided power spectrum
    Pss(2:end,:) = 2*Pss(2:end,:); % double values of single-sided power spectrum except for DC
    F = 0:df:(Fs/2 - df); % frequency vector corresponding to single-sided power spectrum [Hz]
    B = 1/NFFT * sum((windowFn/alpha).^2); % noise power bandwidth of the window function
    
    % crop frequency range to match sensitivity file and compute time vector
    minF = min(S(:,1)); % find lowest frequency in sensitivity curve [Hz]
    maxF = min([max(S(:,1)) max(F)]); % find highest common frequency between sensitivity curve and FFT frequency vector [Hz]
    [~,minFidx] = min(abs(F - minF)); % find index of FFT frequency vector corresponding to minF
    [~,maxFidx(1)] = min(abs(S(:,1) - maxF)); % find index of sensitivity file frequency corresponding to maxF
    [~,maxFidx(2)] = min(abs(F - maxF)); % find index of FFT frequency vector corresponding to maxF
    S = S(1:maxFidx(1),:); % crop sensitivity to maximum common frequency limit
    Pss = Pss(minFidx:maxFidx(2),:); % crop power spectrum to maximum common frequency limit
    F = F(minFidx:maxFidx(2)); % crop frequency vector to maximum common frequency limit
    dt = (1-overlapFactor)/df; % time step [s]
    T = 0:dt:((size(Pss,2)*dt) - dt); % time vector [s]
    
    % build SpectData structure
    SpectData.time = T;
    SpectData.frequency = F;

    % compute time average of power spectrum, with window length of Tavg
    Pss_Tavg = [];
    numAvg = floor(size(Pss,2)*dt/Tavg); % number of averaging windows
    if numAvg >= 1
        for i = 1:numAvg
            startIdx = round(1 + (i-1)*Tavg/dt); % start time index of window
            endIdx = round(startIdx + Tavg/dt - 1); % end time index of window
            Pss_Tavg(:,i) = mean(Pss(:,startIdx:endIdx),2); % time average of power spectrum (not in dB units at this point which is why we can just take the mean without un-logging)
        end
    else
        Pss_Tavg = NaN(size(Pss,1),1);
    end
    
    % compute power spectrum (PS) or power spectral density (PSD) depending on specified specType and assign to SpectData structure
    if strcmp(specType,'power')
        PS = 10*log10(Pss) - S(:,2); % power spectrum in decibels [dB re 1µPa^2]
        PSavg = 10*log10(mean(10.^(PS/10),2)); % time-averaged power spectrum over full file duration [dB re 1µPa^2] 
        PSavgT = 10*log10(Pss_Tavg) - S(:,2); % time-averaged power spectrum in windows of duration Tavg [dB re 1µPa^2] 
        SpectData.P = PS;
        SpectData.Pavg = PSavg;
        SpectData.PavgT = PSavgT;
        SpectData.specType = 'Power Spectrum'; % record that P is 'power spectrum'
        SpectData.specUnits = 'dB re 1µPa^2'; % units of P
        SpectData.Tavg = Tavg; % time average duration [s]
    elseif strcmp(specType,'psd')
        % note that the factor 1/(B*df) corrects for the noise power bandwidth of the window function B, and standardizes the results to the levels which would result from a 1 Hz bin spacing by dividing by the frequency bin spacing df      
        PSD = 10*log10(1/(B*df) * Pss) - S(:,2); % power spectral density in decibels [dB re 1µPa^2/Hz]
        PSDavg = 10*log10(mean(10.^(PSD/10),2)); % time-averaged power spectral density over full file duration [dB re 1µPa^2/Hz] 
        PSDavgT = 10*log10(1/(B*df) * Pss_Tavg) - S(:,2); % time-averaged power spectrum in windows of duration Tavg [dB re 1µPa^2/Hz] 
        SpectData.P = PSD;
        SpectData.Pavg = PSDavg;
        SpectData.PavgT = PSDavgT;
        SpectData.specType = 'Power Spectral Density'; % record that P is 'power spectral density'
        SpectData.specUnits = 'dB re 1µPa^2/Hz'; % units of P
        SpectData.Tavg = Tavg; % time average duration [s]
    end
    SpectData.processingComment = ['Hann window with ' num2str(100*overlapFactor) '% overlap, ' num2str(NFFT) ' samples per FFT window, ' ...
            'temporal resolution: ' num2str(dt) ' s, spectral resolution: ' num2str(df) ' Hz. Sample rate: ' num2str(Fs) ' Hz. ' ...
            'Calibration is applied to each frequency bin individually.'];
end