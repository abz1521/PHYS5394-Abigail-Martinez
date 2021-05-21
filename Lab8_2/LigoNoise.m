% Ligo Noise Simulation
% Read Ligo Sensitivity PSD
LIGOSen = load('iLIGOSensitivity.txt', '-ascii');

% Pre-filtering
lowCutoff = 50; % Hz
highCuttoff = 700; % Hz

% Low Frequency Cutoff
lowSn = LIGOSen(42,2); 

for i = 1:42
    LIGOSen(i,2) = lowSn;
end

% High Frequency Cutoff
highSn = LIGOSen(71,2);

for i = 71:length(LIGOSen(:,1))
    LIGOSen(i,2) = highSn;
end

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
freqVals  = dataLIGO(:,1);
PSDvals   = dataLIGO(:,2);
filtrOrdr = 600; 

b = fir2(fltrOrdr,freqVals/(samplFreq/2),sqrt(PSDvals));

noiseLIGOReal = sqrt(samplFreq)*fftfilt(b,noise);

%% Estimate the PSD using pwelch

[pxx, f] = pwelch(noiseLIGOReal, samplFreq, [], [], sampFreq);

figure;
loglog(f,pxx);
title('Estimated PSD of LIGO Noise Realization');
xlabel('Frequency (kHz)');
ylabel('PSD');

figure;
plot(noiseLIGOReal);
