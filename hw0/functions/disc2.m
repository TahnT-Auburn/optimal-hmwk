%% Probem 6 - Compare Discrete/Continuous Sims
function [sys_obsvd, sys_contrd, sys_obsv_contrd, sys_compd, sys_comp_cld] = ...
disc2( sys_obsv, sys_contr, sys_obsv_contr, sys_comp, sys_comp_cl, T)

% discrete observer
sys_obsvd = c2d(sys_obsv, T, 'zoh');

% discrete controller
sys_contrd = c2d(sys_contr, T, 'zoh');

% discrete observer/controller
sys_obsv_contrd = c2d(sys_obsv_contr, T, 'zoh');

% discrete compensator
sys_compd = c2d(sys_comp, T, 'zoh');

% discrete closed loop compensator
sys_comp_cld = c2d(sys_comp_cl, T, 'zoh');

end