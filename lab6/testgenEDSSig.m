% Signal parameters
snr = 10;
ta = 0.3;
f = 20;
tau = 0.2;
phi = 0;
l = 0.6;

maxFreq = 20;

samplFreq = 1024; %hz
samplInterval = 1/samplFreq;
x = 0:samplInterval:1.0; % 1 second in 1024 resolution

sigVec = genEDSSig(x,snr,ta,f,tau,phi,l);

%% Plot the spectrogram
winLen = 0.05;
ovrlp  = 0.04;

% Convert to integer number of samples
winLenSampl = floor(winLen*samplFreq);
ovrlpSampl  = floor(ovrlp*samplFreq);
[S,F,T]     = spectrogram(sigVec, winLenSampl, ovrlpSampl,[],samplFreq);

% Plot the spectrogram
figure;
imagesc(T,F,abs(S)); axis xy;
title('Exponentially damped Sinusoid:')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
