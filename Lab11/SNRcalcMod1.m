%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

%%
% Folder containing signal and noise generation codes has been added to
% repository

%%
% This is the target SNR for the LR
snr = 10;

%%
% Data generation parameters
nSamples = 2048;
sampFreq = 1024;
x = (0:(nSamples-1))/sampFreq;


%%
% Generate the Exponentially Damped Sinusoidal Signal 
ta = 0.3;
f = 20;
tau = 0.2;
phi = 0;
l = 0.6;
% Amplitude value does not matter as it will be changed in the normalization
A = 1; 
sigVec = genEDSSig(x,snr,ta,f,tau,phi,l);

%%
% We will use the noise PSD used in colGaussNoiseDemo.m but add a constant
% to remove the parts that are zero. (Exercise: Prove that if the noise PSD
% is zero at some frequencies but the signal added to the noise is not,
% then one can create a detection statistic with infinite SNR.)
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

%%
% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies. 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,psdPosFreq);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
%Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdPosFreq);
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,psdPosFreq);
end
%%
% Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
legend('H_0','H_1');
title(['Estimated SNR = ',num2str(estSNR)]);


%% Generate and plot the periodogram of signal and data

% FFT of the data
fftData = fft(dataVec);
% Take positive fourier frequencies
fftData = fftData(1:kNyq);
% FFT of signal (positive frequencies)
fftSig = fft(sigVec);
fftSig = fftSig(1:kNyq);

% Plot periodogram of signal
figure;
plot(posFreq, abs(fftData));
hold on;
plot(posFreq, abs(fftSig));
title('Periodogram of Signal and Data');
xlabel('Frequency (Hz)');
ylabel('Periodogram');

%% Generate and plot the spectrogram

% Set the window length and overlap
winLen = 0.05; %s
ovrlp  = 0.04; %s

winLenSampl = floor(winLen*sampFreq);
ovrlpSampl  = floor(ovrlp*sampFreq);
[S,F,T]     = spectrogram(dataVec, winLenSampl, ovrlpSampl, [], sampFreq);

% Plot spectrogram
figure;
imagesc(T, F, abs(S)); axis xy;
title('Spectrogram of Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure;
plot(x,dataVec);
hold on;
plot(x,sigVec);
xlabel('Time (sec)');
ylabel('Data');