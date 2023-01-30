%% Problem 4 - Equivalent Compensator
function sys_comp = comp(K, L, sys)

% system matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% compensator matrices
A_comp = A - B*K - L*C;
B_comp = L;
C_comp = K;
D_comp = D;

% compensator system
sys_comp = ss(A_comp, B_comp, C_comp, D_comp);

end