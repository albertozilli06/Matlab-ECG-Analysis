% clear all windows
close all;
clear all;

% load ECG data
ECG = load('ECG.mat', 'EKG1', 'EKG2', 'EKG3', 'EKG4', 'EKG5');

EK1 = ECG.EKG1;
EK2 = ECG.EKG2;
EK3 = ECG.EKG3;
EK4 = ECG.EKG4;
EK5 = ECG.EKG5;

% Figure 1: Unfiltered ECGs
figure
p1 = subplot(3, 2, 1);
plot([1:length(EK1)]./250, EK1); % transforming the axes into time
title('Original EKG1')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p2 = subplot(3, 2, 2);
plot([1:length(EK2)]./250, EK2);
title('Original EKG2')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p3 = subplot(3, 2, 3);
plot([1:length(EK3)]./360, EK3);
title('Original EKG3')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p4 = subplot(3, 2, 4);
plot([1:length(EK4)]./360, EK4);
title('Original EKG4')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p5 = subplot(3, 2, 5);
plot([1:length(EK5)]./250, EK5);
title('Original EKG5')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

axis([p1,p2,p3,p4,p5],[0 10 -inf inf]) % displaying time from 0-10s with adjustable y axis
axis 'auto y'

%% Filtering

% cutting off the signal as the T waves become bigger than the R peaks
% at x=295s, where the T wave will be as high as the R peak
EK5 = EK5(1:73750); % 295s * 250 samples/s

D = [4 40]/(250);
F = [4 19]/360;
[y,z] = butter(3,D); % bandpass filter
[t,u] = butter(3,F); % bandpass filter

EK1ST = filtfilt(y, z, EK1);
EK2ST = filtfilt(y, z, EK2);
EK3ST = filtfilt(t, u, EK3);
EK4ST = filtfilt(t, u, EK4);
EK5ST = filtfilt(y, z, EK5);

% Figure 2: Filtered ECGs
figure
p11 = subplot(3, 2, 1);
plot([1:length(EK1ST)]./250, EK1ST);
title('Bandpass Filtered EKG1')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p21 = subplot(3, 2, 2);
plot([1:length(EK2ST)]./250, EK2ST);
title('Bandpass Filtered EKG2')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p31 = subplot(3, 2, 3);
plot([1:length(EK3ST)]./360, EK3ST);
title('Bandpass Filtered EKG3')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p41 = subplot(3, 2, 4);
plot([1:length(EK4ST)]./360, EK4ST);
title('Bandpass Filtered EKG4')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

p51 = subplot(3, 2, 5);
plot([1:length(EK5ST)]./250, EK5ST);
title('Bandpass Filtered EKG5')
xlabel('Time [s]')
ylabel('Amplitude [mV]')

axis([p11,p21,p31,p41,p51],[0 10 -inf inf])
axis 'auto y'


%% Peak detection

% R peak detection

% only detects one peak every 0.2s
[y1ST, x1ST] = findpeaks(EK1ST, 'MINPEAKDISTANCE', round(0.2 * 250));
[y2ST, x2ST] = findpeaks(EK2ST, 'MINPEAKDISTANCE', round(0.2 * 250));
[y3ST, x3ST] = findpeaks(EK3ST, 'MINPEAKDISTANCE', round(0.2 * 360));
[y4ST, x4ST] = findpeaks(EK4ST, 'MINPEAKDISTANCE', round(0.2 * 360));
[y5ST, x5ST] = findpeaks(EK5ST, 'MINPEAKDISTANCE', round(0.2 * 250));


% Calculating average as threshold above which all other R peaks should be

peak_average_1 = 0;
for i=1:length(y1ST)
    peak_average_1 = peak_average_1 + y1ST(i);
end
peak_average_1 = peak_average_1/length(y1ST);

peak_average_2 = 0;
for i=1:length(y2ST)
    peak_average_2 = peak_average_2 + y2ST(i);
end
peak_average_2 = peak_average_2/length(y2ST);

peak_average_3 = 0;
for i=1:length(y3ST)
    peak_average_3 = peak_average_3 + y3ST(i);
end
peak_average_3 = peak_average_3/length(y3ST);

peak_average_4 = 0;
for i=1:length(y4ST)
    peak_average_4 = peak_average_4 + y4ST(i);
end
peak_average_4 = peak_average_4/length(y4ST);

peak_average_5 = 0;
for i=1:length(y5ST)
    peak_average_5 = peak_average_5 + y5ST(i);
end
peak_average_5 = peak_average_5/length(y5ST);
% has to be modified as the T wave grows as big as the R peak with time
peak_average_5 = peak_average_5*1.65;


% Detecting R peaks above threshold

xR1ST=[];
yR1ST=[];
inc = 1;
for i = 1:length(x1ST)
    tmp = EK1ST(x1ST(i));
  if i > 0 && tmp > peak_average_1
      yR1ST(inc) = EK1ST(x1ST(i));
      xR1ST(inc) = x1ST(i);
     inc = inc + 1;
  end
end

xR2ST=[];
yR2ST=[];
inc = 1;
for i = 1:length(x2ST)
    tmp = EK2ST(x2ST(i));
  if i > 0 && tmp > peak_average_2
      yR2ST(inc) = EK2ST(x2ST(i));
      xR2ST(inc) = x2ST(i);
     inc = inc + 1;
  end
end

xR3ST=[];
yR3ST=[];
inc = 1;
for i = 1:length(x3ST)
    tmp = EK3ST(x3ST(i));
  if i > 0 && tmp > peak_average_3
      yR3ST(inc) = EK3ST(x3ST(i));
      xR3ST(inc) = x3ST(i);
     inc = inc + 1;
  end
end

xR4ST=[];
yR4ST=[];
inc = 1;
for i = 1:length(x4ST)
    tmp = EK4ST(x4ST(i));
  if i > 0 && tmp > peak_average_4
      yR4ST(inc) = EK4ST(x4ST(i));
      xR4ST(inc) = x4ST(i);
     inc = inc + 1;
  end
end

xR5ST=[];
yR5ST=[];
inc = 1;
for i = 1:length(x5ST)
    tmp = EK5ST(x5ST(i));
  if i > 0 && tmp > peak_average_5
      yR5ST(inc) = EK5ST(x5ST(i));
      xR5ST(inc) = x5ST(i);
     inc = inc + 1;
  end
end


% S peak detection

xS1 = [];
yS1 = [];
yPeaks_S1 = [];
xPeaks_S1 = [];
EK1ST_inv = -EK1ST;
inc = 1;

% only to second last R peak, as the last R peak might have been detected 
% but no S peak is following any more
for i = 1:length(xR1ST)-1 
    % start and end of interval after R peak
    start_interval = xR1ST(i) + 1; % next sample after R peak
    end_interval = xR1ST(i) + 20; % 20 = 0.08s*250 samples/s (approx)
    
    % Look for local minimum
    [yPeaks_S1, xPeaks_S1] = findpeaks(EK1ST_inv(start_interval:end_interval));
    if (length(xPeaks_S1) > 0)
        xS1(inc) = xPeaks_S1(1)+start_interval;
        yS1(inc) = -yPeaks_S1(1);
        inc = inc+1;
    end
end

xS2 = [];
yS2 = [];
yPeaks_S2 = [];
xPeaks_S2 = [];
EK2ST_inv = -EK2ST;
inc = 1;

for i = 1:length(xR2ST)-1
    % start and end of interval after R peak
    start_interval = xR2ST(i) + 1; % next sample after R peak
    end_interval = xR2ST(i) + 20; % 20 = 0.08s*250 samples/s (approx)
    
    % Look for local minimum
    [yPeaks_S2, xPeaks_S2] = findpeaks(EK2ST_inv(start_interval:end_interval));
    if (length(xPeaks_S2) > 0)
        xS2(inc) = xPeaks_S2(1)+start_interval;
        yS2(inc) = -yPeaks_S2(1);
        inc = inc+1;
    end
end

xS3 = [];
yS3 = [];
yPeaks_S3 = [];
xPeaks_S3 = [];
EK3ST_inv = -EK3ST;
inc = 1;

for i = 1:length(xR3ST)-1
    % start and end of interval after R peak
    start_interval = xR3ST(i) + 1; % next sample after R peak
    end_interval = xR3ST(i) + 30; % 30 = 0.08s*360 samples/s (approx)
    
    % Look for local minimum
    [yPeaks_S3, xPeaks_S3] = findpeaks(EK3ST_inv(start_interval:end_interval));
    if (length(xPeaks_S3) > 0)
        xS3(inc) = xPeaks_S3(1)+start_interval;
        yS3(inc) = -yPeaks_S3(1);
        inc = inc+1;
    end
end

xS4 = [];
yS4 = [];
yPeaks_S4 = [];
xPeaks_S4 = [];
EK4ST_inv = -EK4ST;
inc = 1;

for i = 1:length(xR4ST)-1
    % start and end of interval after R peak
    start_interval = xR4ST(i) + 1; % next sample after R peak
    end_interval = xR4ST(i) + 40; % 30 = 0.08s*360 samples/s (approx)
    
    % Look for local minimum
    [yPeaks_S4, xPeaks_S4] = findpeaks(EK4ST_inv(start_interval:end_interval));
    if (length(xPeaks_S4) > 0)
        xS4(inc) = xPeaks_S4(1)+start_interval;
        yS4(inc) = -yPeaks_S4(1);
        inc = inc+1;
    end
end

xS5 = [];
yS5 = [];
yPeaks_S5 = [];
xPeaks_S5 = [];
EK5ST_inv = -EK5ST;
inc = 1;

for i = 1:length(xR5ST)-1
    % start and end of interval after R peak
    start_interval = xR5ST(i) + 1; % next sample after R peak
    end_interval = xR5ST(i) + 20; % 20 = 0.08s*250 samples/s (approx)
    
    % Look for local minimum
    [yPeaks_S5, xPeaks_S5] = findpeaks(EK5ST_inv(start_interval:end_interval));
    if (length(xPeaks_S5) > 0)
        xS5(inc) = xPeaks_S5(1)+start_interval;
        yS5(inc) = -yPeaks_S5(1);
        inc = inc+1;
    end
end


% T wave detection

xT1 = [];
yT1 = [];
yPeaks_T1 = [];
xPeaks_T1 = [];
inc = 1;

for i = 1:length(xR1ST)-1
    % start and end of interval after R peak
    start_interval = xR1ST(i) + 40; % 40 = 0.16s*250 samples/s
    end_interval = xR1ST(i) + 75; % 75 = 0.3s*250 samples/s
    
    % Look for local maximum
    [yPeaks_T1, xPeaks_T1] = findpeaks(EK1ST(start_interval:end_interval));
    if (length(xPeaks_T1) > 0)
        xT1(inc) = xPeaks_T1(1)+start_interval;
        yT1(inc) = yPeaks_T1(1);
        inc = inc+1;
    end
end

xT2 = [];
yT2 = [];
yPeaks_T2 = [];
xPeaks_T2 = [];
inc = 1;

for i = 1:length(xR2ST)-1
    % start and end of interval after R peak
    start_interval = xR2ST(i) + 40; % 40 = 0.16s*250 samples/s
    end_interval = xR2ST(i) + 75; % 75 = 0.3s*250 samples/s
    
    % Look for local maximum
    [yPeaks_T2, xPeaks_T2] = findpeaks(EK2ST(start_interval:end_interval));
    if (length(xPeaks_T2) > 0)
        xT2(inc) = xPeaks_T2(1)+start_interval;
        yT2(inc) = yPeaks_T2(1);
        inc = inc+1;
    end
end

xT3 = [];
yT3 = [];
yPeaks_T3 = [];
xPeaks_T3 = [];
inc = 1;

for i = 1:length(xR3ST)-1
    % start and end of interval after R peak
    start_interval = xR3ST(i) + 110; % trial and error
    end_interval = xR3ST(i) + 150; % trial and error
    
    % Look for local maximum
    [yPeaks_T3, xPeaks_T3] = findpeaks(EK3ST(start_interval:end_interval));
    if (length(xPeaks_T3) > 0)
        xT3(inc) = xPeaks_T3(1)+start_interval;
        yT3(inc) = yPeaks_T3(1);
        inc = inc+1;
    end
end

xT4 = [];
yT4 = [];
yPeaks_T4 = [];
xPeaks_T4 = [];
inc = 1;

for i = 1:length(xR4ST)-1
    % start and end of interval after R peak
    start_interval = xR4ST(i) + 85; % trial and error
    end_interval = xR4ST(i) + 150; % trial and error
    
    % Look for local maximum
    [yPeaks_T4, xPeaks_T4] = findpeaks(EK4ST(start_interval:end_interval));
    if (length(xPeaks_T4) > 0)
        xT4(inc) = xPeaks_T4(1)+start_interval;
        yT4(inc) = yPeaks_T4(1);
        inc = inc+1;
    end
end

xT5 = [];
yT5 = [];
yPeaks_T5 = [];
xPeaks_T5 = [];
inc = 1;

for i = 1:length(xR5ST)-1
    % start and end of interval after R peak
    start_interval = xR5ST(i) + 41; % 40 = 0.16s*250 samples/s
    end_interval = xR5ST(i) + 75; % 75 = 0.3s*250 samples/s
    
    % Look for local maximum
    [yPeaks_T5, xPeaks_T5] = findpeaks(EK5ST(start_interval:end_interval));
    if (length(xPeaks_T5) > 0)
        xT5(inc) = xPeaks_T5(1)+start_interval;
        yT5(inc) = yPeaks_T5(1);
        inc = inc+1;
    end
end


% Figure 3: ST filtered signal and R peaks/S peaks/T waves
figure
p12 = subplot(3, 2, 1);
plot([1:length(EK1ST)]./250, EK1ST);
title('R, S and T in EKG1');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
stem(xS1./250, yS1, 'g'); % displays S peaks
stem(xT1./250, yT1, 'o'); % displays T wave max
stem(xR1ST./250, yR1ST, 'k'); % displays R peaks
% plot([0 xR1ST(end)], [peak_average_1 peak_average_1]); % plots peak average of EK1ST
hold off

p22 = subplot(3, 2, 2);
plot([1:length(EK2ST)]./250, EK2ST);
title('R, S and T in EKG2');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
stem(xS2./250, yS2, 'g'); % displays S peaks
stem(xT2./250, yT2, 'o'); % displays T wave max
stem(xR2ST./250, yR2ST, 'k'); % displays R peaks
% plot([0 xR2ST(end)], [peak_average_2 peak_average_2]); % plots peak average
hold off

p32 = subplot(3, 2, 3);
plot([1:length(EK3ST)]./360, EK3ST);
title('R, S and T in EKG3');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
stem(xS3./360, yS3, 'g'); % displays S peaks
stem(xT3./360, yT3, 'o'); % displays T wave max
stem(xR3ST./360, yR3ST, 'k'); % displays R peaks
% plot([0 xR3ST(end)], [peak_average_3 peak_average_3]); % plots peak average
hold off

p42 = subplot(3, 2, 4);
plot([1:length(EK4ST)]./360, EK4ST);
title('R, S and T in EKG4');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
stem(xS4./360, yS4, 'g'); % displays S peaks
stem(xT4./360, yT4, 'o'); % displays T wave max
stem(xR4ST./360, yR4ST, 'k'); % displays R peaks
% plot([0 xR4ST(end)], [peak_average_4 peak_average_4]); % plots peak average
hold off

p52 = subplot(3, 2, 5);
plot([1:length(EK5ST)]./250, EK5ST);
title('R, S and T in EKG5');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
stem(xS5./250, yS5, 'g'); % displays S peaks
stem(xT5./250, yT5, 'o'); % displays T wave max
stem(xR5ST./250, yR5ST, 'k'); % displays R peaks
% plot([0 xR5ST(end)], [peak_average_5 peak_average_5]); % plots peak average
hold off

axis([p12,p22,p32,p42,p52],[0 10 -inf inf])
axis 'auto y'


%% ST Shift

% calculating the middle value between S and T

xST1_middle = [];
yST1_middle = [];
yST1_average = [];

for i = 1:length(xS1)

    % middle value
    xST1(i) = xT1(i) - xS1(i); % difference between S and T
    xST1_middle(i) = xS1(i) + round(xST1(i)/2);
    yST1_middle(i) = EK1ST(xST1_middle(i));   

    % average
    yST1_avrg = 0;
    inc = 0;
    for j = xS1(i):xT1(i)
            yST1_avrg = yST1_avrg + EK1ST(j);
            inc = inc + 1;
    end
    yST1_avrg = yST1_avrg/inc;
    
    yST1_average(i) = yST1_avrg;

end


xST2_middle = [];
yST2_middle = [];
yST2_average = [];

for i = 1:length(xS2)

    % middle value
    xST2(i) = xT2(i) - xS2(i); % difference between S and T
    xST2_middle(i) = xS2(i) + round(xST2(i)/2);
    yST2_middle(i) = EK2ST(xST2_middle(i));

    % average
    yST2_avrg = 0;
    inc = 0;
    for j = xS2(i):xT2(i)
            yST2_avrg = yST2_avrg + EK2ST(j);
            inc = inc + 1;
    end
    yST2_avrg = yST2_avrg/inc;
    
    yST2_average(i) = yST2_avrg;

end


xST3_middle = [];
yST3_middle = [];
yST3_average = [];

for i = 1:length(xS3)

    % middle value
    xST3(i) = xT3(i) - xS3(i); % difference between S and T
    xST3_middle(i) = xS3(i) + round(xST3(i)/2);
    yST3_middle(i) = EK3ST(xST3_middle(i));

    % average
    yST3_avrg = 0;
    inc = 0;
    for j = xS3(i):xT3(i)
            yST3_avrg = yST3_avrg + EK3ST(j);
            inc = inc + 1;
    end
    yST3_avrg = yST3_avrg/inc;
    
    yST3_average(i) = yST3_avrg;
    

end


xST4_middle = [];
yST4_middle = [];
yST4_average = [];

for i = 1:length(xS4)

    % middle value
    xST4(i) = xT4(i) - xS4(i); % difference between S and T
    xST4_middle(i) = xS4(i) + round(xST4(i)/2);
    yST4_middle(i) = EK4ST(xST4_middle(i));

    % average
    yST4_avrg = 0;
    inc = 0;
    for j = xS4(i):xT4(i)
            yST4_avrg = yST4_avrg + EK4ST(j);
            inc = inc + 1;
    end
    yST4_avrg = yST4_avrg/inc;
    
    yST4_average(i) = yST4_avrg;
    

end


xST5_middle = [];
yST5_middle = [];
yST5_average = [];

for i = 1:length(xS5)

    % middle value
    xST5(i) = xT5(i) - xS5(i); % difference between S and T
    xST5_middle(i) = xS5(i) + round(xST5(i)/2);
    yST5_middle(i) = EK5ST(xST5_middle(i));

    % average
    yST5_avrg = 0;
    inc = 0;
    for j = xS5(i):xT5(i)
            yST5_avrg = yST5_avrg + EK5ST(j);
            inc = inc + 1;
    end
    yST5_avrg = yST5_avrg/inc;
    
    yST5_average(i) = yST5_avrg;
    
end


% Figure 4: ST middle value over time
figure
subplot(3, 2, 1);
plot([1:length(xST1_middle)], yST1_middle, 'r');
title('ST Shift over time in EKG1');
xlabel('Time');
ylabel('Amplitude [mV]');
hold on
plot([1:length(xST1_middle)], yST1_average, 'g');
hold off
% legend('Middle values between S and T', 'Average between S and T');

subplot(3, 2, 2);
plot([1:length(xST2_middle)], yST2_middle, 'r');
title('ST Shift over time in EKG2');
xlabel('Time');
ylabel('Amplitude [mV]');
hold on
plot([1:length(xST2_middle)], yST2_average, 'g');
hold off
% legend('Middle values between S and T', 'Average between S and T');

subplot(3, 2, 3);
plot([1:length(xST3_middle)], yST3_middle, 'r');
title('ST Shift over time in EKG3');
xlabel('Time');
ylabel('Amplitude [mV]');
hold on
plot([1:length(xST3_middle)], yST3_average, 'g');
hold off
% legend('Middle values between S and T', 'Average between S and T');

subplot(3, 2, 4);
plot([1:length(xST4_middle)], yST4_middle, 'r');
title('ST Shift over time in EKG4');
xlabel('Time');
ylabel('Amplitude [mV]');
hold on
plot([1:length(xST4_middle)], yST4_average, 'g');
hold off
% legend('Middle values between S and T', 'Average between S and T');

subplot(3, 2, 5);
plot([1:length(xST5_middle)], yST5_middle, 'r');
title('ST Shift over time in EKG5');
xlabel('Time');
ylabel('Amplitude [mV]');
hold on
plot([1:length(xST5_middle)], yST5_average, 'g');
hold off
% legend('Middle values between S and T', 'Average between S and T');


% Figure 5: Middle value and average over EKGs
figure
p14 = subplot(3, 2, 1);
plot([1:length(EK1ST)]./250, EK1ST);
title('ST Shift in EKG1');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
plot(xST1_middle./250, yST1_middle, 'r');
plot(xST1_middle./250, yST1_average, 'g');
% legend('EKG1 filtered', 'Middle values between S and T', 'Average between S and T');
% stem(xST1_middle./250, yST1_middle, 'r');
% stem(xST1_middle./250, yST1_average, 'g');
hold off

p24 = subplot(3, 2, 2);
plot([1:length(EK2ST)]./250, EK2ST);
title('ST Shift in EKG2');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
plot(xST2_middle./250, yST2_middle, 'r');
plot(xST2_middle./250, yST2_average, 'g');
% legend('EKG2 filtered', 'Middle values between S and T', 'Average between S and T');
% stem(xST2_middle./250, yST2_middle, 'r');
% stem(xST2_middle./250, yST2_average, 'g');
hold off

p34 = subplot(3, 2, 3);
plot([1:length(EK3ST)]./360, EK3ST);
title('ST Shift in EKG3');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
plot(xST3_middle./360, yST3_middle, 'r');
plot(xST3_middle./360, yST3_average, 'g');
% legend('EKG3 filtered', 'Middle values between S and T', 'Average between S and T');
% stem(xST3_middle./360, yST3_middle, 'r');
% stem(xST3_middle./360, yST3_average, 'g');
hold off

p44 = subplot(3, 2, 4);
plot([1:length(EK4ST)]./360, EK4ST);
title('ST Shift in EKG4');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
plot(xST4_middle./360, yST4_middle, 'r');
plot(xST4_middle./360, yST4_average, 'g');
% legend('EKG4 filtered', 'Middle values between S and T', 'Average between S and T');
% stem(xST4_middle./360, yST4_middle, 'r');
% stem(xST4_middle./360, yST4_average, 'g');
hold off

p54 = subplot(3, 2, 5);
plot([1:length(EK5ST)]./250, EK5ST);
title('ST Shift in EKG5');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
hold on
plot(xST5_middle./250, yST5_middle, 'r');
plot(xST5_middle./250, yST5_average, 'g');
% legend('EKG5 filtered', 'Middle values between S and T', 'Average between S and T');
% stem(xST5_middle./250, yST5_middle, 'r');
% stem(xST5_middle./250, yST5_average, 'g');
hold off

axis([p14,p24,p34,p44,p54],[0 10 -inf inf])
axis 'auto y'