%% Optimal Homework 3 - Problem 1
% Tahn Thawainin
clc
clear variables
close all

%% Simulation Specs

% simulation time
t_sim = 100;

% sampling rate
dt = 0.1;

%% Closed Loop Continuos System

A = [0, 1;
     -1, -1.4];

B = [0; 1];

C = [1, 0];

D = 0;

%% Noise Paramters

% process noise
Q = 4;

% measurement noise
R = 1;

%% Discrete Approximation

% discrete state-transition via matrix exponential
Ad = expm(A*dt);

% discrete input matrix
Bd = dt.*B;

% discrete process noise
Qd = B*Q*B'*dt;

% alternative method
s = [-A, B*Q*B';
     zeros(2,2), A'];

c = expm(s*dt);

% alt discrete state-transition
Ad2 = [c(3,3:4); c(4,3:4)]';

% alt discrete process noise 
Qd2 = [c(3,3:4); c(4,3:4)]'*[c(1,3:4); c(2,3:4)];

%% Simulate System

% time
t = 0:dt:t_sim;

% noise input
w = sqrt(Q)*randn(length(t),1);

% simulation observation
C_sim = [1, 0;
         0, 1];

% continuos system
sys = ss(A,B,C_sim,D);

% cont sim
Yc = lsim(sys,w,t);

% disc sim
Yd = dlsim(Ad,Bd,C_sim,D,w);

% manual simulation
% initialize
x_dd = 0;
x_d = 0;
x = 0;

for i = 1:length(t)
    
    % acceleration
    x_dd = A(2,1)*x + A(2,2)*x_d + w(i);
    
    % integrate for velocity and position
    x_d = x_d + x_dd*dt;
    x = x + x_d*dt;

    % siphon variables
    accel(i) = x_dd;
    vel(i) = x_d;
    pos(i) = x;
end

%% Measurements

% position measurements
y_pos = pos + sqrt(Q)*randn(1,length(t));

%% Noise Paramters

% process noise (ideal ~ 3)
Q_KF = 0.01;

% measurement noise (ideal ~ 5.75)
R_KF = 1;

%% Iterative Discrete KF

% initialize
x_init = [0; 0];
P_init = [.5, 0;
          0, .5];

x = x_init;
P = P_init;

for i = 1:length(t)

    % time update
    x = Ad*x;
    
    P = Ad*P*Ad' + B*Q_KF*B'*dt;
    
    % track priori covariance
    P_pri(:,:,i) = P;

    % measurement update
    K = P*C'/(C*P*C' + R_KF);

    P = (eye(2) - K*C)*P;
    
    % track posteriori covariance
    P_post(:,:,i) = P;

    x = x + K*(y_pos(i) - C*x);
    
    % estimate errors
    pos_err(i) = x(1) - pos(i);
    vel_err(i) = x(2) - vel(i);

    % track variables
    pos_est(i) = x(1);
    vel_est(i) = x(2);
    kalmanGain(:,:,i) = K;
 

end

%% Error Analysis

% standard deviations
err_p_std = std(pos_err);
err_v_std = std(vel_err);

% norm
N = sqrt(err_v_std^2 + err_p_std^2);lting 


%% Steady-State Analysis

for i = 2:length(t)
    
    if mod(kalmanGain(:,:,i),kalmanGain(:,:,i-1)) == 0

        % steady state kalman gain
        Kss = kalmanGain(:,:,i);
    else
        continue
    end

    if mod(P_pri(:,:,i), P_pri(:,:,i-1)) == 0

        % steady state priori covariance
        P_pri_ss = P_pri(:,:,i);
    else
        continue
    end

    if mod(P_post(:,:,i), P_post(:,:,i-1)) == 0

        % steady state priori covariance
        P_post_ss = P_post(:,:,i);
    else
        continue
    end

end

% closed loop estimator
A_KF = Ad - Kss*C;

% estimator eigenvalues
eig_KF = eig(A_KF);

%% Interface
% 0 - false
% 1 - true

% display simulated system
disp_sim = 1;
% display estimator specs
disp_specs = 1;
% display steady state info
disp_ss = 1;
% display solution plot
disp_sol = 1;

% simulation plots
if disp_sim == 1

    figure
    subplot(2,1,1)
    plot(t, pos, LineWidth = 1.5)
    title('Simulated Position and Velocity')
    ylabel('Position')
    subplot(2,1,2)
    plot(t,vel, LineWidth = 1.5)
    ylabel('Velocity')
    xlabel('Time [s]')
    set(gcf, 'color', 'w')

elseif disp_sim == 0
end

% estimator specs
if disp_specs == 1

    disp('Process Noise')
    disp(Q_KF)
    disp('Measurement Noise')
    disp(R_KF)

elseif disp_specs == 0
end

% steady state info
if disp_ss == 1

    disp('Steady State Kalman Gain:')
    disp(Kss)
    disp('Steady State Prior Covariance')
    disp(P_pri_ss)
    disp('Steady State Posterirori Covariance')
    disp(P_post_ss)
    disp('Estimator Eigenvalues')
    disp(eig_KF)

    % kalman gains plots
    figure
    hold on
    plot(t,kalmanGain(1,:), LineWidth=2)
    plot(t,kalmanGain(2,:), LineWidth=2)
    hold off
    legend('K_1', 'K_2')
    title('Kalman Gains')
    xlabel('Time [s]')
    set(gcf, 'color', 'w')

elseif disp_ss == 0
end

% KF/error info
if disp_sol == 1

    disp('Error STDs:')
    disp('Position Error STD')
    disp(err_p_std)
    disp('Velocity Error STD')
    disp(err_v_std)
    disp('Error Norm:')
    disp(N)

    figure
    subplot(3,2,[1,2])
    hold on
    plot(t, pos, '--',DisplayName='Truth',LineWidth=1.5)
    plot(t, pos_est, DisplayName='KF', LineWidth=1.5)
    plot(t, y_pos, '.', DisplayName='Measurements')
    hold off
    ylabel('Position')
    xlabel('Time [s]')
    title('Kalman Filter Position and Velocity Estimates and Errors')
    legend
    set(gcf, 'color', 'w') 
    subplot(3,2,[3,4])
    hold on
    plot(t, vel, '--',DisplayName='Truth',LineWidth=1.5)
    plot(t, vel_est, DisplayName='KF',LineWidth=1.5)
    ylabel('Velocity')
    set(gcf, 'color', 'w')
    
    % error plots    
    subplot(3,2,5)
    plot(t, pos_err, LineWidth=1.5)
    ylabel('Position Error')
    set(gcf, 'color', 'w')    
    subplot(3,2,6)
    plot(t, vel_err, LineWidth=1.5)
    ylabel('Velocity Error')
    set(gcf, 'color', 'w')

elseif disp_sol == 0
end


    















