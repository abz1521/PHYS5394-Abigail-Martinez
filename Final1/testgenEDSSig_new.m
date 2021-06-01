P = struct('tau', 0.2, 'ta', 0.3,'freq', 20,'phi', 0, 'length', 0.6);
           
snr = 10; % amplitude

% Generate time vector of 1 second
samplFreq = 1024; %hz
samplInterval = 1/samplFreq;
dataX = 0:samplInterval:1.0; % 1 second in 1024 resolution

sigVec = genEDSSig_new(dataX,snr,P);

% Plot the signal
figure;
plot(dataX, sigVec, 'Marker', '.', 'MarkerSize', 6);
title('Exponentially Damped Sinusoidal Signal Using Struct');
xlabel('Seconds (s)');
ylabel('Amplitude');

% Generate Periodogram
nSampl  = length(dataX); % # of samples
Ldata= dataX(end)-dataX(1); % seconds of data
% kNyq is the frequency sample corresponding to the Nyquist frequency
kNyq = floor(nSampl/2)+1;
posFreq = (0:(kNyq-1))*(1/Ldata); % vector of positive Fourier frequencies

fftSig = fft(sigVec);
fftSig = fftSig(1:kNyq);

% Plot Periodogram
figure;
plot(posFreq, abs(fftSig));
title('Periodogram of Exponentially Damped Sinusoidal Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');