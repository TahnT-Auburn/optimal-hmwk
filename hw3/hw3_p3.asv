%% Optimal Homework 3 - Problem 3
% Tahn Thawainin

clc
clear variables
% close all

%% Load Data

% data file
data = dlmread('hw3_3.txt', ' ', 2, 0);

% fix data
% NOTE: 
% col1: time, col2: East, col3: North, col4: Psi, col5: gyro, col6: radar 
data = [data(:,4), data(:,6), data(:,9), data(:,12), data(:,15), data(:,18)];

%% Measurements

% time
t = data(:,1);

% smapling rate
dt = 1/5;

% East (m)
E = data(:,2);

% North (m)
N = data(:,3);

% Psi (rad)
Psi = data(:,4);

% gyro (yaw rate - rad/s)
Psi_d = data(:,5);

% radar (velocity - m/s)
v = data(:,6);

%% Noise Parameters

% process noise
Qd = [1, 0, 0, 0, 0, 0, 0;
      0, 1, 0, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0, 0;
      0, 0, 0, 1, 0, 0, 0;
      0, 0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 0, 0.0001, 0;
      0, 0, 0, 0, 0, 0, 0.0001;];

Q_rls = 0;

% measurement noise
Rd =      [1, 0, 0, 0, 0;
           0, 1, 0, 0, 0;
           0, 0, 1, 0, 0;
           0, 0, 0, 1, 0;
           0, 0, 0, 0, 1];

%% Kalman Filter/Recursive LS

% filter type
% 0 - RLS
% 1 - KF
filterType = 1;

% outage condition
% 0 - false
% 1 - true
outage = 1;

% initialize
x_init = [E(1); N(1); Psi(1); v(1); Psi_d(1); 0; 0];

P_init = [1, 0, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0, 0;
          0, 0, 1, 0, 0, 0, 0;
          0, 0, 0, 1, 0, 0, 0;
          0, 0, 0, 0, 1, 0, 0;
          0, 0, 0, 0, 0, 1, 0;
          0, 0, 0, 0, 0, 0, 1;];

x = x_init;
P = P_init;

for i = 1:length(t)
    
    % state transition matrix
    Ad = [1, 0, 0, sin(Psi(i))*dt, 0, -sin(Psi(i))*dt, 0;
          0, 1, 0, cos(Psi(i))*dt, 0, -cos(Psi(i))*dt, 0;
          0, 0, 1, 0, dt, 0, -dt;
          0, 0, 0, 1, 0, 0, 0;
          0, 0, 0, 0, 1, 0, 0;
          0, 0, 0, 0, 0, 1, 0;
          0, 0, 0, 0, 0, 0, 1];
    
%     % input matrix
%     Bd = [sin(Psi(i))*dt, 0;
%           cos(Psi(i))*dt, 0;
%           0, dt;
%           0, 0;
%           0, 0];
%     
%     % inputs
%     u = [v(i);
%          Psi_d(i)];
    
    % observation matrix
    C =  [1, 0, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0, 0;
          0, 0, 1, 0, 0, 0, 0;
          0, 0, 0, 1, 0, 1, 0;
          0, 0, 0, 0, 1, 0, 1];
    
    % measurements
    y = [E(i);
         N(i);
         Psi(i);
         v(i);
         Psi_d(i)];

    % RLS filter
    if filterType == 0

    % RLS update
    K = P*C'/(C*P*C' + Rd);

    P = (eye(7) - K*C)*P;
    
    x = x + K*(y - C*x);
    
    % KF filter - no outage
    elseif filterType == 1 && outage == 0
    
    % time update
    x = Ad*x; %+ Bd*u;
    
    P = Ad*P*Ad' + Qd;

    % measurement update
    K = P*C'/(C*P*C' + Rd);

    P = (eye(7) - K*C)*P;
    
    x = x + K*(y - C*x);
    
    % KF filter - outage
    elseif filterType == 1 && outage == 1
    
    if t(i) < 260
    % time update
    x = Ad*x; %+ Bd*u;
    
    P = Ad*P*Ad' + Qd;

    % measurement update
    K = P*C'/(C*P*C' + Rd);

    P = (eye(7) - K*C)*P;
    
    x = x + K*(y - C*x);

    % outage at last 40 seconds
    elseif t(i) >= 260
    % time update
    x = Ad*x; %+ Bd*u;
    
    P = Ad*P*Ad' + Qd;
    end

    end

    % siphon variables
    E_est(i) = x(1);
    N_est(i) = x(2);
    Psi_est(i) = x(3);
    v_est(i) = x(4);
    Psi_d_est(i) = x(5);
    Bias_r(i) = x(6);
    Bias_g(i) = x(7);

    kalmanGain(:,:,i) = K;
    covar(:,:,i) = P;
end

%% Least Squares

% measurement batch
y_ls = [E';
        N';
        Psi';
        v';
        Psi_d'];

for k = 1:length(t)

% least squares solution
x_ls(:,k) = pinv(C)*y_ls(:,k);

end
% siphon variables
E_est_ls = x_ls(1,:)';
N_est_ls = x_ls(2,:);
Psi_est_ls = x_ls(3,:);
v_est_ls = x_ls(4,:);
Psi_d_est_ls = x_ls(5,:);
Bias_r_ls = x_ls(6,:);
Bias_g_ls = x_ls(7,:);

%% Means

% radar measurement mean
v_mean = mean(v);

% gyro measurement mean
Psi_d_mean = mean(Psi_d);

% estimated radar velocity mean
v_est_mean = mean(v_est);

% estimated gyro yaw rate mean
Psi_d_est_mean = mean(Psi_d_est);

% LS means
v_est_mean_ls = mean(v_est_ls);
Psi_d_mean_ls = mean(Psi_d_est_ls);

% estimated radar bias mean (KF or RLS and LS)
bias_r_mean = mean(Bias_r);
bias_r_mean_ls = mean(Bias_r_ls);

% estimated gyro bias mean (KF or RLS and LS)
bias_g_mean = mean(Bias_g);
bias_g_mean_ls = mean(Bias_g_ls);

%% Interface
% 0 - false
% 1 - true

% RGB Triplets
darkBlue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
lightBlue = [0.3010 0.7450 0.9330];
red = [0.6350 0.0780 0.1840];

% display KF solution
disp_KF_RLS_sol = 1;
% display LS solution
disp_LSsol = 1;
% display estimator performance
disp_est_perf = 1;

% KF solution info
if disp_KF_RLS_sol == 1 && filterType == 1

    disp('Filter Type: Kalman Filter')

    figure
    subplot(3,1,1)
    hold on
    plot(t, E_est, DisplayName='KF', LineWidth=1.75)
    plot(t, E, '--', DisplayName='Meas', LineWidth=1.75)
    hold off
    if outage == 1
        xlim([250, 300])
    else
    end
    ylabel('meters')
    title('East Estimate')
    legend
    set(gcf, 'color', 'w')    
    subplot(3,1,2)
    hold on
    plot(t, N_est, DisplayName='KF', LineWidth=1.75)
    plot(t, N, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    if outage == 1
        xlim([250, 300])
    else
    end
    ylabel('meters')
    title('North Estimate')
    set(gcf, 'color', 'w')
    subplot(3,1,3)
    hold on
    plot(t, Psi_est, DisplayName='KF', LineWidth=1.75)
    plot(t, Psi, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    if outage == 1
        xlim([250, 300])
    else
    end
    xlabel('Time [s]')
    ylabel('rad')
    title('Yaw Estimate')
    set(gcf, 'color', 'w')
    
    figure
    subplot(2,2,1)
    plot(t, Bias_g, LineWidth=1.75)
    ylabel('Gyro Bias (rad/s)')
    xlabel('Time [s]')
    title('Gyro Bias')
    subplot(2,2,2)
    plot(t, Bias_r, LineWidth=1.75)
    title('Radar Bias')
    ylabel('Radar Bias (m/s)')
    xlabel('Time [s]')
    set(gcf, 'color', 'w')
    
    subplot(2,2,3)
    hold on
    plot(t, Psi_d, LineWidth=1.75, DisplayName='Meas')
    plot(t, Psi_d_est, LineWidth=1.75, DisplayName='KF')
    hold off
    legend
    ylabel('$\dot{\phi}$ (rad/s)', 'interpreter','latex')
    xlabel('Time [s]')
    title('Gyro Meas VS Gyro Estimate')
    set(gcf, 'color', 'w')
    ylim([-0.05, 0.1])
    subplot(2,2,4)
    hold on
    plot(t, v, LineWidth=1.75, DisplayName='Meas')
    plot(t, v_est, LineWidth=1.75, DisplayName='KF')
    hold off
    ylabel('$v$ (m/s)', 'interpreter','latex')
    xlabel('Time[s]')
    title('Radar Meas VS Radar Estimate')
    ylim([2, 2.3])
    set(gcf, 'color', 'w')

% RLS solution
elseif disp_KF_RLS_sol == 1 && filterType == 0
    
    disp('Filter Type: Recursive Least Squares')
    
    figure
    subplot(3,1,1)
    hold on
    plot(t, E_est, DisplayName='RLS', LineWidth=1.75)
    plot(t, E, '--', DisplayName='Meas', LineWidth=1.75)
    hold off
    ylabel('meters')
    title('East Estimate')
    legend
    set(gcf, 'color', 'w')    
    subplot(3,1,2)
    hold on
    plot(t, N_est, DisplayName='RLS', LineWidth=1.75)
    plot(t, N, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    ylabel('meters')
    title('North Estimate')
    set(gcf, 'color', 'w')
    subplot(3,1,3)
    hold on
    plot(t, Psi_est, DisplayName='RLS', LineWidth=1.75)
    plot(t, Psi, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    xlabel('Time [s]')
    ylabel('rad')
    title('Yaw Estimate')
    set(gcf, 'color', 'w')
    
    figure
    subplot(2,2,1)
    plot(t, Bias_g, LineWidth=1.75)
    ylabel('Gyro Bias (rad/s)')
    xlabel('Time [s]')
    title('Gyro Bias')
    subplot(2,2,2)
    plot(t, Bias_r, LineWidth=1.75)
    title('Radar Bias')
    ylabel('Radar Bias (m/s)')
    xlabel('Time [s]')
    set(gcf, 'color', 'w')
    
    subplot(2,2,3)
    hold on
    plot(t, Psi_d, LineWidth=1.75, DisplayName='Meas')
    plot(t, Psi_d_est, LineWidth=1.75, DisplayName='RLS')
    hold off
    legend
    ylabel('$\dot{\phi}$ (rad/s)', 'interpreter','latex')
    xlabel('Time [s]')
    title('Gyro Meas VS Gyro Estimate')
    set(gcf, 'color', 'w')
    ylim([-0.05, 0.1])
    subplot(2,2,4)
    hold on
    plot(t, v, LineWidth=1.75, DisplayName='Meas')
    plot(t, v_est, LineWidth=1.75, DisplayName='RLS')
    hold off
    ylabel('$v$ (m/s)', 'interpreter','latex')
    xlabel('Time[s]')
    title('Radar Meas VS Radar Estimate')
    ylim([2, 2.3])
    set(gcf, 'color', 'w')

elseif disp_KF_RLS_sol == 0
end

% LS solution info
if disp_LSsol == 1
    
    figure
    subplot(3,1,1)
    hold on
    plot(t, E_est_ls, DisplayName='LS', LineWidth=1.75)
    plot(t, E, '--', DisplayName='Meas', LineWidth=1.75)
    hold off
    ylabel('meters')
    title('East Estimate')
    legend
    set(gcf, 'color', 'w')    
    subplot(3,1,2)
    hold on
    plot(t, N_est_ls, DisplayName='LS', LineWidth=1.75)
    plot(t, N, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    ylabel('meters')
    title('North Estimate')
    set(gcf, 'color', 'w')
    subplot(3,1,3)
    hold on
    plot(t, Psi_est_ls, DisplayName='LS', LineWidth=1.75)
    plot(t, Psi, '--', DisplayName='Meas',LineWidth=1.75)
    hold off
    xlabel('Time [s]')
    ylabel('rad')
    title('Yaw Estimate')
    set(gcf, 'color', 'w')
    
    figure
    subplot(2,2,1)
    plot(t, Bias_g_ls, LineWidth=1.75)
    ylabel('Gyro Bias (rad/s)')
    xlabel('Time [s]')
    title('Gyro Bias')
    subplot(2,2,2)
    plot(t, Bias_r_ls, LineWidth=1.75)
    title('Radar Bias')
    ylabel('Radar Bias (m/s)')
    xlabel('Time [s]')
    set(gcf, 'color', 'w')
    
    subplot(2,2,3)
    hold on
    plot(t, Psi_d, LineWidth=1.75, DisplayName='Meas')
    plot(t, Psi_d_est_ls, LineWidth=1.75, DisplayName='LS')
    hold off
    legend
    ylabel('$\dot{\phi}$ (rad/s)', 'interpreter','latex')
    xlabel('Time [s]')
    title('Gyro Meas VS Gyro Estimate')
    set(gcf, 'color', 'w')
    ylim([-0.05, 0.1])
    subplot(2,2,4)
    hold on
    plot(t, v, LineWidth=1.75, DisplayName='Meas')
    plot(t, v_est_ls, LineWidth=1.75, DisplayName='LS')
    hold off
    ylabel('$v$ (m/s)', 'interpreter','latex')
    xlabel('Time[s]')
    title('Radar Meas VS Radar Estimate')
    ylim([2, 2.3])
    set(gcf, 'color', 'w')

elseif disp_LSsol == 0
end

if disp_est_perf == 1
    
    % KF/RLS estimated biases
    disp('Gyro Bias Mean (absolute value)')
    disp(abs(bias_g_mean))
    disp('Radar Bias Mean (absolute value)')
    disp(abs(bias_r_mean))
    
    disp('Gyro Meas/Est mean difference')
    disp(abs(Psi_d_mean - Psi_d_est_mean))
    disp('Radar Meas/Est mean differience')
    disp(abs(v_mean - v_est_mean))

    % LS estimated biases
    disp('LS Gyro Bias Mean (AV)')
    disp(abs(bias_g_mean_ls))
    disp('LS Radar Bias Mean (AV)')
    disp(abs(bias_r_mean_ls))

    disp('Gyro Meas/Est mean difference')
    disp(abs(Psi_d_mean - Psi_d_mean_ls))
    disp('Radar Meas/Est mean differience')
    disp(abs(v_mean - v_est_mean_ls))

elseif disp_est_perf == 0
end












