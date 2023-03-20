%% Optimal Homework 2
% Tahn Thawainin

clear
clc
close all

%% Problem 3 - Simulate a gyroscope and estimate scale and bias factors

% samples
simp = 1000; 

% set scale and bias factors
a = 1.3;
b = 2;

% set time
t = 1/simp:1/simp:1;

% gyro dynamics
w = 2*pi;              % frequency
r = 100*sin(w*t);

% measurement covariance
R = deg2rad(0.3)^2;

% % noise
% n = deg2rad(0.3)*randn(1,length(r));
% 
% % gyro measurements
% g = a*r + b + n;
% 
% % least squares
% 
% % observation matrix
% H = [r', ones(simp,1)];
% %H = [r((i-1)*simp+1:i*simp)', ones(simp,1)];
% 
% % DOP
% DOP = diag(inv(H'*H));
% 
% % predicted state covariance
% P_est = R*DOP;
% 
% % measurement
% y = g';
% %y = g((i-1)*simp+1:i*simp)';
% 
% % state estimate
% x = (H'*H)\H'*y;

% monte carlo
for k = 1:1000        

    % noise
    n = deg2rad(0.3)*randn(1,length(r));

    % gyro measurements
    g = a*r + b + n;

    % least squares

    % observation matrix
    H = [r', ones(simp,1)];
    %H = [r((i-1)*simp+1:i*simp)', ones(simp,1)];
    
    % DOP
    DOP = diag(inv(H'*H));
    
    % predicted state covariance
    P_est = R*DOP;

    % measurement
    y = g';
    %y = g((i-1)*simp+1:i*simp)';
    
    % state estimate
    x = (H'*H)\H'*y;

    % monte carlo results
    x_amc_ls(k,:) = x(1,:);
    x_bmc_ls(k,:) = x(2,:);

    x_amc_errors_ls(k,:) = x_amc_ls(k,:) - a;
    x_bmc_errors_ls(k,:) = x_bmc_ls(k,:) - b;
end

% monte carlo error mean and std
x_acoeff_mean_ls = mean(x_amc_ls);
x_bcoeff_mean_ls = mean(x_bmc_ls);

x_means_lsa = mean(x_amc_errors_ls);
x_means_lsb = mean(x_bmc_errors_ls);

x_std_lsa = std(x_amc_errors_ls);
x_std_lsb = std(x_bmc_errors_ls);

%% recursive least squares

% set rls time
t = 1/simp:1/simp:100;

% gyro dynamics
r_rls = 100*sin(w*t);

% initialize
x_init = [0.5; 0.5];

P_init = [10, 0;
          0, 10];
x = x_init;
P = P_init;

% monte carlo
for kk = 1:100

    % noise
    n = deg2rad(0.3)*randn(1,length(r_rls));

    % gyro measurements
    g = a*r_rls + b + n;

    % initialize
    x_init = [0.5; 0.5];
    
    P_init = [10, 0;
              0, 10];
    x = x_init;
    P = P_init;
    
    % propagate 
    for ii = 1:50
        
        % observation matrix
        H = [r_rls((ii-1)*simp+1:ii*simp)', ones(simp,1)];
        
        % DOP
        DOP_rls(ii,:) = diag(inv(H'*H));
        
        % predicted state covariance
        P_est_rls(ii,:) = R*DOP_rls(ii,:);
        
        % measurement
        y = g((ii-1)*simp+1:ii*simp)';
        
        % pre-fit residual covariance
        S = H*P*H' + R*(eye(simp));
    
        % gain
        L = P*H'/S;
    
        % state estimate
        x = x + L*(y - H*x);

        % update state covariance
        P = (eye(2) - L*H)*P;

        % variables
        x_rls(ii,:) = x;
        
        P_rls(:,:,ii) = P;

    end

    % monte carlo results
    x_amc_rls(:,kk) = x_rls(:,1);
    x_bmc_rls(:,kk) = x_rls(:,2);

    x_amc_errors_rls(:,kk) = x_amc_rls(:,kk) - a;
    x_bmc_errors_rls(:,kk) = x_bmc_rls(:,kk) - b;

end

% monte carlo error mean and std
x_acoeff_mean_rls = mean(x_amc_rls,2);
x_bcoeff_mean_rls = mean(x_bmc_rls,2);

x_acoeff_std_rls = std(x_amc_rls,0,2);
x_bcoeff_std_rls = std(x_bmc_rls,0,2);

x_means_rlsa = mean(x_amc_errors_rls,2);
x_means_rlsb = mean(x_bmc_errors_rls,2);

x_std_rlsa = std(x_amc_errors_rls,0,2);
x_std_rlsb = std(x_bmc_errors_rls,0,2);



%% plotting
% theoretical standard deviation
stda = sqrt(P_rls(1,1,:));
stda = reshape(stda, 1,50);
stdb = sqrt(P_rls(2,2,:));
stdb = reshape(stdb, 1,50);
grayColor = [.7 .7 .7];

% monte carlo results
figure
subplot(2,1,1)
title('Monte Carlo Estimates and 95% Confidence Intervals')
hold on
plot(x_acoeff_mean_rls + 3.*x_acoeff_std_rls, 'r', LineWidth=2)
plot(x_acoeff_mean_rls + 3.*stda', '--k', LineWidth=2)
plot(x_acoeff_mean_rls - 3.*x_acoeff_std_rls, 'r', LineWidth=2)
plot(x_acoeff_mean_rls - 3.*stda', '--k', LineWidth=2)
plot(x_amc_rls ,'color', grayColor)
hold off
ylabel('a Estimates')
legend('Experimental STD', 'Theoretical STD')
subplot(2,1,2)
hold on
plot(x_bcoeff_mean_rls + 3.*x_bcoeff_std_rls, 'r', LineWidth=2)
plot(x_bcoeff_mean_rls + 3.*stdb', '--k', LineWidth=2)
plot(x_bcoeff_mean_rls - 3.*x_bcoeff_std_rls, 'r', LineWidth=2)
plot(x_bcoeff_mean_rls - 3.*stdb', '--k', LineWidth=2)
plot(x_bmc_rls ,'color', grayColor)
ylabel('b Estimates')
xlabel('RLS Runs')
set(gcf,'Color','w')

% error analyses
figure
subplot(2,1,1)
title('Monte Carlo Estimate Errors and 95% Confidence Intervals')
hold on
plot(x_means_rlsa + 3.*x_std_rlsa, 'r', LineWidth=2)
plot(x_means_rlsa + 3.*stda', '--k', LineWidth=2)
plot(x_means_rlsa - 3.*x_std_rlsa, 'r', LineWidth=2)
plot(x_means_rlsa - 3.*stda', '--k', LineWidth=2)
plot(x_amc_errors_rls ,'color', grayColor)
hold off
legend('Experimental STD', 'Theoretical STD')
ylabel('a Estimates Errors')
subplot(2,1,2)
hold on
plot(x_means_rlsb + 3.*x_std_rlsb, 'r', LineWidth=2)
plot(x_means_rlsb + 3.*stdb', '--k', LineWidth=2)
plot(x_means_rlsb - 3.*x_std_rlsb, 'r', LineWidth=2)
plot(x_means_rlsb - 3.*stdb', '--k', LineWidth=2)
plot(x_bmc_errors_rls ,'color', grayColor)
ylabel('b Estimates Errors')
xlabel('RLS Runs')
set(gcf,'Color','w')





