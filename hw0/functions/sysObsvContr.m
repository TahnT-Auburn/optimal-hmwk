%% Observer/Controller Closed Loop System
function sys_obsv_contr = sysObsvContr(sys, sys_contr, sys_obsv, K)

% system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% closed loop controller A matrix 
A_contr = sys_contr.A;

% closed loop observer A matrix
A_obsv = sys_obsv.A;

% closed Loop sys - combined controller and observer
Acl = [A_contr, B*K;
       zeros(size(A)), A_obsv];

Bcl = [B;
       zeros(size(B))];

Ccl = [C, zeros(size(C))];

Dcl = D;

% Closed loop system - combined controller & observer
sys_obsv_contr = ss(Acl,Bcl,Ccl,Dcl);

end 