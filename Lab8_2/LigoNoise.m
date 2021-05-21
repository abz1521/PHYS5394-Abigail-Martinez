% Ligo Noise Simulation
% Read Ligo Sensitivity PSD
LIGOSen = load('iLIGOSensitivity.txt', '-ascii');

% Pre-filtering
lowCutoff = 50; % Hz
highCuttoff = 700; % Hz

% Low Frequency Cutoff
lowSn = LIGOSen(42,2); 

LIGOSen(1:42,2) = lowSn;

% High Frequency Cutoff
highSn = LIGOSen(71,2);

LIGOSen(71:end,2) = highSn;

LIGOSen(2:end+1,:) = LIGOSen;
LIGOSen(1,1) = 0; %Hz
LIGOSen(1,2) = lowSn;

% Sample Frequecy at 6 kHz
samplFreq = 6000;

LIGOSen = LIGOSen(1:85,:); % Cut window to 85
LIGOSen(85,1) = 3000; %Hz

% Generate White Gaussian Noise realization and pass through LIGO PSD filter
nSampl = 20*samplFreq; % 2 (s) * sampFreq
% Gaussian Noise Vector
noise = randn(1,nSampl);

% FIR Filter
freqVals  = LIGOSen(:,1);
PSDvals   = LIGOSen(:,2);
filtrOrdr = 600; 

b = fir2(fltrOrdr,freqVals/(samplFreq/2),sqrt(PSDvals));

noiseLIGOReal = sqrt(samplFreq)*fftfilt(b,noise);

%% Estimate the PSD using pwelch

[pxx, f] = pwelch(noiseLIGOReal, samplFreq, [], [], samplFreq);

figure;
loglog(f,pxx);
title('Estimated PSD of LIGO Noise Realization');
xlabel('Frequency (kHz)');
ylabel('PSD');

figure;
plot(noiseLIGOReal);
