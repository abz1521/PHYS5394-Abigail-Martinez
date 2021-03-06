% Read testData.txt
path = "testData.txt";

data = load(path);

timeVec = data(:,1); % time vector
dataVal = data(:,2); % Data Vector

samplFreq = length(timeVec)/round(timeVec(end)); % Sample Frequency Calculation

% Focuses on the first 5 seconds
colGaus = dataVal(1:5*1024); 

% Estimate PSD for the the first 5 seconds using pwelch
% closGaus gets first 5 seconds, 
%80 is hamming window length
[pxx,f] = pwelch(colGaus, 80,[],[],samplFreq); 

% Design Whitening Filter
b = fir2(400, f / (1024 / 2), 1./sqrt(pxx));
% Whiten TestData
whitenedData = sqrt(samplFreq)*fftfilt(b,dataVal);

% Spectrograms 
winLen = 0.05; % window length
ovrlp = 0.04; % overlap length
 
winLenSamp = floor(winLen*samplFreq);
ovrlpSamp  = floor(ovrlp*samplFreq);

% Generate Spectrogram 
[S,F,T] = spectrogram(dataVal, winLenSamp, ovrlpSamp, [], samplFreq);
% Plot Spectrogram
figure;
imagesc(T, F, abs(S)); axis xy;
title('Spectrogram of Test Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Generate Whitened Spectrogram
[S1,F1,T1] = spectrogram(whitenedData, winLenSamp, ovrlpSamp, [], samplFreq);
% Plot Whitened Spectrogram
figure;
imagesc(T1, F1, abs(S1)); axis xy;
title('Spectrogram of Whittened Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;


figure;
subplot(2,1,1), plot(timeVec, dataVal), title('Test Data');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2), plot(timeVec, whitenedData), title('Whitened Data');
xlabel('Time (s)');
ylabel('Amplitude')

