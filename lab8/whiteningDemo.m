path = "testData.txt";
data1 = load(path);
  
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
title('Spectrogram of testData');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;