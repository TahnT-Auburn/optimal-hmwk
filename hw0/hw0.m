%% Optimal Homework 0
% Tahn Thawainin

clc; clear; close all

%% Model System (Problem 1)

% rotational moment of inertia
J = 10; % kgm^2
% rotational damping
b = 1; % Nms/rad

% system matrices
A = [0,    1;
     0, -b/J];

B = [0; 1/J];

C = [1, 0];

D = 0;

% state space system
sys = ss(A, B, C, D);

% eigenvalues
eigc = eig(A);

%% Call functions

%% Observer function
[O, L, p_obsv, sys_obsv] = observer(50, 0.7, sys);

% check Observability
if rank(O) == length(A)
    disp('System is observable')
else 
    disp('System is NOT observable')
end

% plot Observer step response
figure
step(sys_obsv)
title('Observer Step Response')
set(gcf,"Color",'w')

%% Controller function
[Co, K, p_contr, sys_contr] = controller(10, 0.7, sys);

% Check Controllability
if rank(Co) == length(A)
    disp('System is controllable')
else 
    disp('System is NOT controllable')
end

%% Combined Observer/Controller
sys_obsv_contr = sysObsvContr(sys, sys_contr, sys_obsv, K);

% plot combined observer/controller
figure
step(sys_obsv_contr);
title('Combined Observer and Controller Response')
set(gcf,"Color",'w')

%% Compensator Function
 sys_comp = comp(K, L, sys);

% closed loop compensator transfer function
disp('compensator transfer function')
tf(sys_comp)

% compensator bode plot
figure
bode(sys_comp)
title('Compensator Bode')
grid
set(gcf,"Color",'w')

%% Closed Loop Compensator
 sys_comp_cl = sysCompCL(sys, sys_comp);

% closed loop compensator transfer function
disp('closed loop compensator transfer function')
tf(sys_comp_cl)

% closed loop compensator bode plot
figure
margin(sys_comp_cl)
title('Closed Loop Compensator Bode')
grid
set(gcf,"Color",'w')

% gain and phase margins
[Gm, Pm, Wcg, Wcp] = margin(sys_comp_cl);
disp('closed loop gain/phase margins (cont)')
disp('Gain margin:')
disp(Gm)
disp('Phase margin:')
disp(Pm)

%% Discretized System

% sample rate
T = 1/1000;

[sysd, eigd, p_obsvd, p_contrd, Ld, Kd] = disc(sys, p_obsv, p_contr, T);

% discrete eigenvalues
disp('discrete eigenvalues:')
disp(eigd)

% discrete observer poles
disp('discrete observer poles:')
disp(p_obsvd)

% discrete controller poles
disp('discrete controller poles:')
disp(p_contr)

% discrete L
disp('discrete L:')
disp(Ld);

% discrete K
disp('discrete K:')
disp(Kd);

%% Simulation Comparisions (Problem 6)

% call discretize function
[sys_obsvd, sys_contrd, sys_obsv_contrd, sys_compd, sys_comp_cld] = ...
disc2(sys_obsv, sys_contr, sys_obsv_contr, sys_comp, sys_comp_cl, T);

% continuous vs discretized closed loop observer sys
figure
hold on
step(sys_obsv)
step(sys_obsvd)
hold off
legend('continuous', 'discrete')
set(gcf,"Color",'w')
title('Observer Step Response')

% % continuous vs discretized closed loop controller sys
% figure
% hold on
% step(sys_contr)
% step(sys_contrd)
% hold off
% legend

% continuous vs discretized closed loop observer and controller sys
figure
hold on
step(sys_obsv_contr)
step(sys_obsv_contrd)
hold off
legend('continuous', 'discrete')
set(gcf,"Color",'w')
title('Combined Controller/Observer Step Response')

% continuous vs discretized compensator sys
figure
hold on
step(sys_comp)
step(sys_compd)
hold off
legend('continuous', 'discrete')
set(gcf,"Color",'w')
title('Compensator Step Response')

% continuous vs discretized compensator bode
figure
hold on
margin(sys_comp)
margin(sys_compd)
hold off
grid
legend('continuous', 'discrete')
set(gcf,"Color",'w')
title('Compensator Bode')

% continuous vs discretized closed loop observer sys
figure
hold on
step(sys_comp_cl)
step(sys_comp_cld)
hold off
legend('continuous', 'discrete')
set(gcf,"Color",'w')
title('Closed Loop Compensator Step Response')

% continuous vs discretized closed loop observer bode
figure
hold on
margin(sys_comp_cl)
margin(sys_comp_cld)
hold off
legend('continuous', 'discrete')
grid
set(gcf,"Color",'w')
title('Closed Loop Compensator Bode')

% discrete compensator transfer function
disp('discrete compensator TF:')
tf(sys_compd)

% discrete closed loop compensator transfer function
disp('discrete closed loop compensator TF:')
tf(sys_comp_cld)

