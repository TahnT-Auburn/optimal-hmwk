%% Optimal Homework 1

clc;
clear variables;
close all;

%% Problem 1 - Dice Rolls

% a) six die numbered {1-6}
pX_a = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];    % PDF of a single die roll
pZ_a = pX_a; % initialize
for i = 1:5
    pZ_a = conv(pX_a,pZ_a);   % PDF of six die roll
end

% b) six die numbered {4,5,6,7,8,9}
pZ_b = pZ_a;

% c) six die numbered {1,1,3,3,3,5}
pX_c = [1/3, 0, 1/2, 0, 1/6];    % PDF of a sing die roll
pZ_c = pX_c; % initialize
for i = 1:5
    pZ_c = conv(pX_c,pZ_c);   % PDF of six die roll
end

% d) three die numbered {1-6}, and three die numbered {1,1,3,3,3,5}
pX_d1 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
pX_d2 = [1/3, 0, 1/2, 0, 1/6];
pZ_d = pX_d1; % initialize
for i = 1:5
    if i < 3
        pZ_d = conv(pX_d1,pZ_d);
    else
        pZ_d = conv(pX_d2,pZ_d);
    end
end

% plot time
t_a = 6:36;
t_b = 24:54;
t_c = 6:30;
t_d = 6:33;
% plot PDFs
figure
hold on
stem(t_a, pZ_a*100, DisplayName='PDF_a')
stem(t_b, pZ_b*100, DisplayName='PDF_b')
stem(t_c, pZ_c*100, DisplayName='PDF_c')
stem(t_d, pZ_d*100, DisplayName='PDF_d')
hold off
title('PDFs of Die Rolls')
xlabel('Sums')
ylabel('Probability (%)')
set(gcf,'Color','w')
grid
legend (Location='best')

%% Problem 2 - Joint Random Variables

% joint PDF for two fair dice
Jpdf = (1/36)*ones(6);

% values of x1 and x2
x1 = [1, 2, 3, 4, 5, 6];
x2 = [1, 2, 3, 4, 5, 6];

% a)
% E(x1)
a1 = mean(x1);
% E(x1 - E(x1))
a2 = mean(x1) - mean(x1);
% E(x1^2)
a3 = sum(x1.^2 .*(1/6));
% E((x1 - E(x1))^2) = var(x)
a4 = sum(x1.^2 .*(1/6))  - 2*mean(x1)*mean(x1) + mean(x1)^2;
% E(((x1 - E(x1))(x2 - E(x2)))
a5 = mean(x1)*mean(x2) - mean(x1)*mean(x2) - mean(x2)*mean(x1) + mean(x1)*mean(x2);

% b) covariance matrix for x1 and x2
P = [a4, 0;
     0, a4];

% v1 and v2
v1 = x1;
v2 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

% pdfs of x1 and x2
px1 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
px2 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];

% pdfs of v1 and v2
pZ_v1 = px1;
pZ_v2 = conv(px1, px2);

% c) % joint PDF for variables v1 and v2
for i = 1:length(pZ_v1)
    jpdf(i,:) = pZ_v1(i).*pZ_v2;

    % check sums
    s(i) = sum(jpdf(i,:));
end

% d) find the mean, E(v1 - E(v1)), RMS, and variance of v1
v1_rms = sqrt((1/6)*sum(v1.^2));
v1_var = sum(v1.^2.*pZ_v1) - mean(v1)^2;

% e) find the mean, E(v2 - E(v2)), RMS, and variance of v2
v2_rms = sqrt((1/11)*sum(v2.^2));
v2_var = sum(v2.^2.*pZ_v2) - mean(v2)^2;

%% Problem 4 - Serially-Correlated Random Sequence

% dice values
d1 = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5];
d2 = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5];

% pdf of each die
pd1 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
pd2 = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];

% a)
% pdf of sums
pzd = conv(pd1, pd2);

% possible sums of d1 and d2
Vo = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5];

% b)
% mean
m = pzd*Vo';

% variance
Vo_var = sum(d1.^2*pd1') + sum(d2.^2*pd2');


mc_length = 1000;              % monte carlo length
r = 0.75;                         % r value
Vn = zeros(100,mc_length);     % initialize Vn
Vn_var_mc = zeros(100,mc_length); % initialize Vn variance
% run monte carlo
for i = 1:mc_length

    for k = 1:100

        % define rolls
        roll1 = randi(6);
        roll2 = randi(6);

        % define dice indices
        Vo_mc = d1(roll1) + d2(roll2);
        
        % define Vn
        Vn(k+1, i) = (1-r)*Vn(k,i) + r*Vo_mc;
        
        % variance function
        Vn_var_mc(k+1, i) = (1-r)^2*Vn_var_mc(k,i) + r^2*Vo_var;
    end
end

% generate mean and variance of Vn
for i = 1:100
    Vn_col_mean = mean(Vn(i,:));
    Vn_col_var = var(Vn(i,:));

    % monte carlo variance
    Vn_col_var_mc = mean(Vn_var_mc(i,:));
end
% mean of Vn
Vn_mean_abs = mean(Vn_col_mean);
Vn_var_abs = mean(Vn_col_var);
Vn_var_mc_abs = mean(Vn_col_var_mc);

figure
plot(Vn)
title('Vn, r = -0.1')
xlabel('Rolls')
ylabel('Vn values')
grid
set(gcf,'color','w')

%% Problem 6

% covariance matrix
Px = [2 1;
      1 4];

% eigenvalue and eigenvector
[V, D] = eig(Px);

% transform matrix
A = D^(-1/2)*V;

% circle matrix
% angles for a cirlce
T = 0:360;
for i = 1:length(T)
a1(i) = cosd(T(i));
a2(i) = sind(T(i));
a(:,i) = [a1(i); a2(i)];
end

% distinct c values
a1 = 0.25*a;
a2 = a;
a3 = 1.5*a;

% ellipses
b1 = inv(A)*a1;
b2 = inv(A)*a2;
b3 = inv(A)*a3;

% plot ellipses
figure
hold on
fill(b1(1,:),b1(2,:),'cyan', 'FaceAlpha', 1, DisplayName='c = 0.25')
fill(b2(1,:),b2(2,:),'m', 'FaceAlpha', 0.5, DisplayName='c = 1')
fill(b3(1,:),b3(2,:),'g', 'FaceAlpha', 0.2, DisplayName='c = 1.5')
hold off
title('Likelihood Ellipses')
set(gcf,'Color','w')
grid
legend (Location='best')

%% Problem 7

% standard deviation
sigmaX = 2;

% generate random vector x & y
x = -10:0.01:10;
y = 2.*(x.^2);

% PDF of x & y
for i = 1:length(x)
pdfX(i) = (1/(sqrt(2*pi)*sigmaX))*exp((-1/2) * (x(i)^2/sigmaX^2));
fX_scale(i) = (1/(sqrt(2*pi)*sigmaX));
fX_expo(i) = exp((-1/2) * (x(i)^2/sigmaX^2));
pdfY(i) = (1/(4*sqrt(pi*y(i))*sigmaX))*exp((-1/2) * (y(i)/(2*sigmaX^2)));
fY_scale(i) = (1/(4*sqrt(pi*y(i))*sigmaX));
fY_expo(i) = exp((-1/2) * (y(i)/(2*sigmaX^2)));
end

figure
hold on
plot(x,pdfX*100, DisplayName='PDF X')
plot(y,pdfY*100, DisplayName='PDF Y')
hold off
title('PDFs of x & y')
xlabel('Values')
ylabel('Probability (%)')
xlim([-50, 50])
ylim([0, 100])
grid
legend(Location='best')
set(gcf,'color','w')

