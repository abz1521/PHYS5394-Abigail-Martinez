trials = 10000;
% Parameters
a = -2;
b = 1;
m = 1.5;
s = 2;

for i = 1:trials;
% Generate the vector of trials from uniform PDF
uniformPDF(i) = customrand(a,b); % Between -2,1

% Generate the vector of trials from normal PDF
normalPDF(i) = customrandn(m,s); % Mean = 1.5, Std = 2.0

end

figure;

subplot(2,1,1); % 2 rows, 1 coloumn, 1 of 2 plots
plot(uniformPDF); % Plotting Uniform PDF
hold on;
title('Uniform PDF');
xlabel('Trial');
ylabel('Value');

subplot(2,1,2); % 2 rows, 1 coloumn, 2 of 2 plots
histogram(uniformPDF, 'Normalization', 'PDF'); % Plotting histogram
title('Histogram');
xlabel('Value');
ylabel('Probability');

figure;

subplot(2,1,1); % 2 rows, 1 coloumn, 1 of 2 plots
plot(normalPDF); % Plotting Gauss PDF
hold on;
title('Gauss PDF');
xlabel('Trial');
ylabel('Value');

subplot(2,1,2); % 2 rows, 1 coloumn, 2 of 2 plots
histogram(normalPDF, 'Normalization', 'PDF'); % Plotting histogram for Gauss PDF
title('Histogram');
xlabel('Value');
ylabel('Probability');

