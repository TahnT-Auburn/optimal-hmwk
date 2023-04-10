%% Optimal Homework 3 - Problem 2
% Tahn Thawainin

clc
clear variables
close all

%% Load Data

% data file
data = dlmread('hw3_2.txt');

% time
t = data(:,1);

% measurement
y = data(:,2);

%% Discrete System
Ad = 1;
B = 0;
C = 1;
D = 0;

%% Noise Parameters

% process noise
Q = 0.01;

% measurement noise
R = 1;

%% Kalman Filter

% initialize
x_init = 0;
P_init = 1;

x = x_init;
P = P_init;
  
for i = 1:length(t)

    % time update
    x = Ad*x;
    
    P = Ad*P*Ad' + Q;

    % measurement update
    K = P*C'/(C*P*C' + R);

    P = (eye(1) - K*C)*P;

    x = x + K*(y(i) - C*x);
    
    % siphon variables
    bias_est(i) = x;
    kalmanGain(i) = K;
 
end

%% Steady-State Analysis

for i = 2:length(t)
    
    if mod(kalmanGain(i),kalmanGain(i-1)) < 0.0001

        % steady state kalman gain
        Kss = kalmanGain(i);
    else
        continue

    end

end

%% Low Pass Filter

% first order low-pass filter 
num = sqrt(Q);
den = [1, -1 + sqrt(Q)];

% filtered solution
yf = filter(num,den,y,y(1));

%% Interface
% 0 - false
% 1 - true

% display estimator specs
disp_specs = 1;
% display solution plot
disp_sol = 0;
% display steady state
disp_ss  = 0;
% display low-pass filter info
disp_lpfilt= 1;

% estimator specs
if disp_specs == 1

    disp('Process Noise')
    disp(Q)
    disp('Measurement Noise')
    disp(R)

elseif disp_specs == 0
end

% KF solution info
if disp_sol == 1

    figure
    hold on
    plot(t, y, '.')
    plot(t, bias_est, LineWidth=1.75)
    hold off
    ylabel('Bias Estimate')
    xlabel('Time [s]')
    title('Kalman Filter Measurement Bias Estimates Using Various Q values')
%     legend('Measurements', 'Q = 0.001')
    set(gcf, 'color', 'w') 

elseif disp_sol == 0
end

% steady state info
if disp_ss == 1  && Q > 0

    disp('Steady State Kalman Gain:')
    disp(Kss)

    % steady state kalman gains plots
    figure
    hold on
    plot(t, kalmanGain, LineWidth=2)
    hold off
    title('Kalman Gain')
    xlabel('Time [s]')
    set(gcf, 'color', 'w') 

elseif disp_ss == 0
end

if disp_lpfilt == 1 

    figure;
    hold on
    ylabel('Bias')
    xlabel('Time [s]')
    title('Low Pass Filter vs Kalman Filter')   
    plot(t, y, '.', DisplayName='Measurements')
    plot(t, bias_est, DisplayName='KF', LineWidth=1.75)
    plot(t, yf, '--', color = 'm', DisplayName='LPF', LineWidth=1.75)
    legend
    set(gcf, 'color', 'w') 

elseif disp_lpfilt
end

