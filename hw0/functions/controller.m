%% Problem 3 - Design a controller
function [Co, K, p_contr, sys_contr] = controller(fn, zeta, sys)

% system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Controllability matrix
Co = ctrb(A,B);

% Error dynamics
wn = fn*2*pi;
sigma = wn*zeta;
wd = wn*sqrt(1-zeta^2); % Damped frequency

% Observer poles - using second order TF
p_contr = [-sigma + wd*i, -sigma - wd*i];

% design controller
K = place(A,B,p_contr);

% closed loop controller matrix
A_contr = A - B*K;

% closed loop sys - controller
sys_contr = ss(A_contr, B, C, D);

end
