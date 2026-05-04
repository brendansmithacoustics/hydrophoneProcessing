%% Example script to compute absolute sound level measurements from Ocean Networks Canada hydrophone data
% Date/Author:  2026-04-23, Brendan Smith
clear all; close all; clc

% ------------------------------- INPUTS-----------------------------------
% directory containing audio files for analysis and sensitivity .txt file
fdir = '.\hydrophoneProcessing\data\';

% path to audio file
audioPath = [fdir '\ICLISTENHF6327_20260423T000000.000Z.flac'];

% path to calibration .txt file
calPath = [fdir '\ICLISTENHF6327_20230824T000000.000Z-hydrophonePreDeploymentCalibration.txt'];

df = 10; % target frequency bin spacing [Hz], sensitivity file will be interpolated over frequency if necessary
overlapFactor = 0.5; % ratio of overlap between adjacent time windows, typically 0 or 0.5
specType = 'psd'; % specify either 'power' for power spectrum [dB re 1µPa^2] or 'psd' for power spectral densitiy [dB re 1µPa^2/Hz]
Tavg = 60; % averaging time [s]
% -------------------------------------------------------------------------

% read audio file
[Y,Fs] = audioread(audioPath); % read audio file scaled to +/- 1. Do not use 'native' argument. The conversion to bits is done manually below.
aui = audioinfo(audioPath); % extract audio information (sample rate, bit depth, etc.)
bitDepth = aui.BitsPerSample; % audio bit depth [bits]
Y = Y * 2^(bitDepth-1); % convert Y to counts [counts]

% read sensitivity (calibration) file
fileID = fopen(calPath,'r');
Sraw = textscan(fileID,'%f %f','Delimiter',',','HeaderLines', 3); % read both columns as floats, and skip header lines (note older sensitivity .txt files only had one header line)
fileStatus = fclose(fileID);
S(:,1) = Sraw{1,1}; % assign frequency to first column of matrix
S(:,2) = Sraw{1,2}; % assign sensitivity to second column of matrix

% compute absolute sound level spectra from ONC hydrophone audio file
[SpectData] = calibrateONChydrophoneSpectra(Y,Fs,S,df,overlapFactor,specType,Tavg);

%% Plot spectrogram and time averaged spectrum
close all

fSize = 16;
fName = 'Helvetica';
lw = 2;

fh = figure;
fh.Position = [100 100 1200 500];
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile
pcolor(SpectData.time,SpectData.frequency,SpectData.P)
shading flat
xlabel('Time [s]')
ylabel('Frequency [Hz]')
colormap turbo
xlim([0 300])
ylim([10 25600])
clim([40 140])
ch = colorbar;
if strcmp(SpectData.specType,'Power Spectral Density')
    ch.Label.String = 'Power spectral density [dB re 1 μPa^{2}/Hz]';
elseif strcmp(SpectData.specType,'Power Spectrum')
    ch.Label.String = 'Power spectrum [dB re [dB re 1 μPa^{2}]';
end
set(gca,'Layer','top')
set(gca,'LineWidth',1,'TickLength',[0.019 0.019]);
set(gca,'YScale','log')
set(gca,'FontSize',fSize)
set(gca,'FontName',fName)

nexttile
semilogx(SpectData.frequency,SpectData.Pavg,'LineWidth',lw)
grid on
xlim([10 25600])
xlabel('Frequency [Hz]')
if strcmp(SpectData.specType,'Power Spectral Density')
    ylabel('Power spectral density [dB re 1 μPa^{2}/Hz]');
elseif strcmp(SpectData.specType,'Power Spectrum')
    ylabel('Power spectrum [dB re [dB re 1 μPa^{2}]');
end
set(gca,'LineWidth',1);
set(gca,'FontSize',fSize)
set(gca,'FontName',fName)