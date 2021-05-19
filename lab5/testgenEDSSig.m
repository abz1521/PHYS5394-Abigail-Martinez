% Parameters

snr = 10;
ta = 0.3;
f = 20;
tau = 0.2;
phi = 0;
l = 0.6;

%Sample Intervals

sampleFreq = 5*f;
sampleInterval=1/sampleFreq;
x=0:sampleInterval:1.0;

%Run Function

sigVec = genEDSSig(x,snr,ta,f,tau,phi,l);

%Plot Function
%plot(x,sigVec);
%title("EDSS Nyquist Sampling Frequency (1/2)");

%% Plot periodogram
%Length of data
dataLen = x(end)-x(1);
% DFT sample corresponding to Nyquist frequency
nSampl = length(x);
kNyq = floor(nSampl/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Plot the periodogram
figure;
plot(posFreq, abs(fftSig));
title('Periodogram of Exponentially Damped Sinosoidal Wave:')

