path = "testData.txt";

timeVec = data(:,1);
dataVals  = data(:,2);

samplFreq = length(timeVec)/round(timeVec(end));

colGaus = dataVals(1:5*1024);
%  
[pxx,f] = pwelch(colGaus, 80,[],[],samplFreq);
%    
b = fir2(400, f / (1024 / 2), 1./sqrt(pxx));
whitenedData = sqrt(samplFreq)*fftfilt(b,dataVals);
%  
winLen = 0.05;
ovrlp = 0.04;
% 
winLenSamp = floor(winLen*sampFreq);
ovrlpSamp  = floor(ovrlp*sampFreq);
% 
[S,F,T] = spectrogram(dataVals, winLenSamp, ovrlpSamp, [], sampFreq);
% 
figure;
imagesc(T, F, abs(S)); axis xy;
title('Spectrogram of Test Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

[S1,F1,T1] = spectrogram(whitenedData, winLenSamp, ovrlpSamp, [], sampFreq);
figure;
imagesc(T1, F1, abs(S1)); axis xy;
title('Spectrogram of Whittened Data');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

figure;
subplot(2,1,1), plot(timeVec, dataVals), title('Test Data');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2), plot(timeVec, whitenedData), title('Whitened Data');
xlabel('Time (s)');
ylabel('Amplitude')

