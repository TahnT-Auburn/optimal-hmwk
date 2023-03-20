%% Optimal Homework 2
% Tahn Thawainin

clc
clear
% close all

%% Problem 4 - Least Squares for System ID

% transfer function
numd = 0.25*[1 -0.8];   % numerator
dend = [1 -1.9 0.95];   % denominator

% random input
u = randn(1000,1);

% simulate discrete sys
y = dlsim(numd, dend, u);

%%
for jj = 1:4

% use distinct sigmas
if jj == 1
    sigma = 0.01;
elseif jj == 2
    sigma = 0.1;
elseif jj == 3
    sigma = 0.5;
elseif jj == 4
    sigma = 1;
else
end

% covariance
R = sigma^2;

% noise
n = sigma*randn(1000,1);

% output w/ noise
Y = y + n;

% generate observation matrix
for i = 3:length(Y)

    H(i,:) = [-Y(i-1), -Y(i-2), u(i-1), u(i-2)];

end
H = H(3:end,:);

% for i = 1:length(Y)-2
%     
%     H(i,:) = [-Y(i+1), -Y(i), u(i+1), u(i)];
% 
% end

% least squares
Y_ls = Y(3:end);
X1 = (H'*H)\(H'*Y_ls);

% simulated transfer function
num_sim = [X1(3,1) X1(4,1)];
den_sim = [1 X1(1,1) X1(2,1)];


% calculate SNR
SNR(jj,:) = snr(y,n);

% track measurement
Y_track(:,jj) = Y;

% track estimates
X_track(:,jj) = X1;

% track sim transfer functions num and den
num_sim_track(jj,:) = num_sim;
den_sim_track(jj,:) = den_sim;
end

% transfer functions
G = tf(numd, dend);
G_sim1 = tf(num_sim_track(1,:), den_sim_track(1,:));
G_sim2 = tf(num_sim_track(2,:), den_sim_track(2,:));
G_sim3 = tf(num_sim_track(3,:), den_sim_track(3,:));
G_sim4 = tf(num_sim_track(4,:), den_sim_track(4,:));


%% Monte Carlo Sim
for j = 1:4

% use distinct sigmas
if j == 1
    sigma = 0.01;
elseif j == 2
    sigma = 0.1;
elseif j == 3
    sigma = 0.5;
elseif j == 4
    sigma = 1;
else
end

R = sigma^2;

% monte carlo
for k = 1:100 
    
    % output w/ noise
    Y = y + sigma*randn(1000,1);
    
    % generate observation matrix
    for i = 3:length(Y)
    
        H(i,:) = [-Y(i-1), -Y(i-2), u(i-1), u(i-2)];
    
    end
    H = H(3:end,:);
    
    % dilution of percision
    DOP = diag(inv(H'*H));
    
    % estimated covariance
    P_est = R*DOP;

    % least squares
    Y_ls = Y(3:end);
    X = (H'*H)\(H'*Y_ls);
    
    % track estimates
    X_mc_a1(k,:) = X(1,1);
    X_mc_a2(k,:) = X(2,1);
    X_mc_b1(k,:) = X(3,1);
    X_mc_b2(k,:) = X(4,1);

    % track DOP and covariance estimate
    DOP_mc(k,:) = DOP;
    P_est_mc(k,:) = P_est;

end

% track estimates
X_a1_track(:,j) = X_mc_a1;
X_a2_track(:,j) = X_mc_a2;
X_b1_track(:,j) = X_mc_b1;
X_b2_track(:,j) = X_mc_b2;


end

% mean and std of estimates (distinct noise parameters)
mean_a11 = mean(X_a1_track(:,1));
mean_a12 = mean(X_a1_track(:,2));
mean_a13 = mean(X_a1_track(:,3));
mean_a14 = mean(X_a1_track(:,4));

mean_a21 = mean(X_a2_track(:,1));
mean_a22 = mean(X_a2_track(:,2));
mean_a23 = mean(X_a2_track(:,3));
mean_a24 = mean(X_a2_track(:,4));

mean_b11 = mean(X_b1_track(:,1));
mean_b12 = mean(X_b1_track(:,2));
mean_b13 = mean(X_b1_track(:,3));
mean_b14 = mean(X_b1_track(:,4));

mean_b21 = mean(X_b2_track(:,1));
mean_b22 = mean(X_b2_track(:,2));
mean_b23 = mean(X_b2_track(:,3));
mean_b24 = mean(X_b2_track(:,4));

std_a11 = std(X_a1_track(:,1));
std_a12 = std(X_a1_track(:,2));
std_a13 = std(X_a1_track(:,3));
std_a14 = std(X_a1_track(:,4));

std_a21 = std(X_a2_track(:,1));
std_a22 = std(X_a2_track(:,2));
std_a23 = std(X_a2_track(:,3));
std_a24 = std(X_a2_track(:,4));

std_b11 = std(X_b1_track(:,1));
std_b12 = std(X_b1_track(:,2));
std_b13 = std(X_b1_track(:,3));
std_b14 = std(X_b1_track(:,4));

std_b21 = std(X_b2_track(:,1));
std_b22 = std(X_b2_track(:,2));
std_b23 = std(X_b2_track(:,3));
std_b24 = std(X_b2_track(:,4));

% mean of theoretical covariance
std_a1_t = sqrt(mean(P_est_mc(:,1)));
std_a2_t = sqrt(mean(P_est_mc(:,2)));
std_b1_t = sqrt(mean(P_est_mc(:,3)));
std_b2_t = sqrt(mean(P_est_mc(:,4)));

%% Plotting

% plot systems
figure
subplot(4,1,1)
title ('Ideal vs Measured Systems with Distinct Noise Parameters')
hold on
plot(y)
plot(Y_track(:,1))
hold off
legend('Ideal', '\sigma = 0.01')
subplot(4,1,2)
hold on
plot(y)
plot(Y_track(:,2))
hold off
legend('Ideal', '\sigma = 0.1')
subplot(4,1,3)
hold on
plot(y)
plot(Y_track(:,3))
hold off
legend('Ideal', '\sigma = 0.5')
subplot(4,1,4)
hold on
plot(y)
plot(Y_track(:,4))
hold off
legend('Ideal', '\sigma = 1')
set(gcf, 'color', 'w')

% plot bode
figure
hold on
margin(G)
margin(G_sim1)
margin(G_sim2)
margin(G_sim3)
margin(G_sim4)
hold off
title('Bode Plots at Distinct Noise Parameters')
legend('Ideal','\sigma = 0.01','\sigma = 0.1', '\sigma = 0.5', '\sigma = 1')
grid
set(gcf, 'color', 'w')






