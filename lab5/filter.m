% Parameters
a1 = 10; %SNR
f1 = 100; % Frequency
phi1 = 0; % Initial Phase

a2 = 5;
f2 = 200;
phi2 = pi/6;

a3 = 2.5;
f3 = 300;
phi3 = pi/4;

% Generate time vector
samplFreq = 1024; 
maxFreq = 0.5*samplFreq; 
nSampl = 2*samplFreq;
timeVec = (0:(nSampl-1))/samplFreq;

% Declaration of Sinusoids
s1 = a1*sin(2*pi*f1*timeVec+phi1);
s2 = a2*sin(2*pi*f2*timeVec+phi2);
s3 = a3*sin(2*pi*f3*timeVec+phi3);

% Sum of Sinusoids
s = s1 + s2 + s3;

%% Filter Designs
% Design low pass filter
filtOrdr = 30;
LowPass = fir1(filtOrdr, 150/maxFreq, 'low');
% Apply lowpass filter
filtsig1 = fftfilt(LowPass,s);

% Filter 2 Design (bandpass)
BandPass = fir1(filtOrdr, [150,250]/maxFreq, 'bandpass');
% Apply Bandpass filter 
filtsig2 = fftfilt(BandPass, s);
filtsig2 = fftfilt(BandPass, filtsig2);

% Filter 3 Design
highPass = fir1(filtOrdr, 250/maxFreq, 'high');
% Apply highpass filter
filtsig3 = fftfilt(highPass, s);

% DFT samplefor Nyquist Frequency
kNyq = floor(nSampl/2)+1;
% Positive Fourier Frequency
posFreq = (0:(kNyq-1))/(nSampl/samplFreq);

% FFT of unfiltered signal
fftSig = fft(s);
fftSig = fftSig(1:kNyq);

% FFT of filtered signal 1
fftSig1 = fft(filtsig1);
fftSig1 = fftSig1(1:kNyq);

% FFT of filtered signal 2
fftSig2 = fft(filtsig2);
fftSig2 = fftSig2(1:kNyq);

% FFT of filtered signal 3
fftSig3 = fft(filtsig3);
fftSig3 = fftSig3(1:kNyq);

% Periodogram Plot 
fig = figure;

subplot(2,2,1); % 2 rows, 2 coloumns, 4 Plots 
hold on;
plot(posFreq,abs(fftSig));
title('Input signal');

subplot(2,2,2); 
hold on;
plot(posFreq,abs(fftSig1));
title('Lowpass');

subplot(2,2,3);
hold on;
plot(posFreq,abs(fftSig2));
title('Bandpass');

subplot(2,2,4);
hold on;
plot(posFreq,abs(fftSig3));
title('Highpass');

han = axes(fig, 'visible', 'off'); % Common axis labels
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'Magnitude');
xlabel(han, 'Frequency (Hz)');

%% Plot the signals
% Filtered
fig2 = figure;

subplot(3,1,1); % 3 row, 1 column, 3 plots
hold on;
plot(timeVec, s);
plot(timeVec, filtsig1);
title('Lowpass cutoff 150 Hz (s1)');

subplot(3,1,2);
hold on;
plot(timeVec, s);
plot(timeVec,filtsig2);
title('Bandpass cutoff 150-250 Hz (s2)');

subplot(3,1,3);
hold on;
plot(timeVec, s);
plot(timeVec, filtsig3);
title('Highpass cutoff 250 Hz (s3)');

han = axes(fig2, 'visible', 'off'); % Common axis label
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'Magnitude');
xlabel(han, 'Time');



