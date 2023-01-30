%% Problem 2 - Design an observer
function [O, L, p_obsv, sys_obsv] = observer(fn, zeta, sys)

% system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% observability matrix
O = obsv(A,C);

% error dynamics
wn = fn*2*pi;
sigma = wn*zeta;
wd = wn*sqrt(1-zeta^2); % damped frequency

% observer poles - using second order TF
p_obsv = [-sigma + wd*i, -sigma - wd*i];

% design observer
L = place(A',C',p_obsv)';

% closed Loop observer matrix
A_obsv = A - L*C;

% closed loop sys - observer
sys_obsv = ss(A_obsv, B, C, D);

end 