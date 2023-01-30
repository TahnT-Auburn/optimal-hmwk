%% Closed Loop Compensator
function sys_comp_cl = sysCompCL(sys, sys_comp)

% system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% compensator matrices
A_comp = sys_comp.A;
B_comp = sys_comp.B;
C_comp = sys_comp.C;
D_comp = sys_comp.D;

% closed loop compensator
Acl = [A, B*C_comp;
        -B_comp*C, A_comp];

Bcl = [zeros(size(B)); B_comp];

Ccl = [C, zeros(size(C))];

Dcl = D;

% closed loop system - compensator
sys_comp_cl = ss(Acl,Bcl,Ccl,Dcl);

end