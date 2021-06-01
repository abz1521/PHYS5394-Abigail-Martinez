% Generate time series
samplFreq = 1024; % Hz
% dataLen = 4.0; % seconds
samplInterval = 1/samplFreq;
dataX = 0:samplInterval:1.0; % time series vector

% Generate the struct of EDSS signal parameters
P = struct('tau', 0.2, 'ta', 0.3,'freq', 20,'phi', 0, 'length', 0.6);

% Function handle to call genSGSig_new
H = @(snr) genEDSSig_new(dataX, snr, P);

%% Plot EDSS signal at 10, 12, and 15 SNR
S = [10 12 15];

legVec = [];
figure;
% Loop through snrs, generate signal and add to plot
for snr = S
    legVec{end+1} = ['snr=',num2str(snr)];
    sigVec = H(snr);
    % Plot the signal
    plot(dataX, sigVec);
    hold on;
end
% Display title, axis labels, and legend
title('Exponentially Damped Sinusoidal Signal at Various SNRs');
xlabel('Time (s)');
ylabel('Amplitude');

hold off;
legend(legVec);
