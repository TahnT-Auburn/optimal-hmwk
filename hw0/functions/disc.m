%% Problem 5 - Discretize System
function [sysd, eigd, p_obsvd, p_contrd, Ld, Kd] = disc(sys, p_obsv, p_contr, T)

% discrete system
sysd = c2d(sys, T, 'zoh');

% system matrices
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

% discrete eigenvalues
eigd = eig(sysd);

% discrete observer poles
p_obsvd = [exp(p_obsv(1)*T), exp(p_obsv(2)*T)];

% discrete observer poles
p_contrd = [exp(p_contr(1)*T), exp(p_contr(2)*T)];

% design discrete observer
Ld = place(Ad',Cd',p_obsvd)';

% design discrete controller
Kd = place(Ad, Bd, p_contrd);

end