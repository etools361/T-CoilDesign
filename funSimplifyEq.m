%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-09-10(yyyy-mm-dd)
% 简化公式
%--------------------------------------------------------------------------
% 简化公式
% 定义符号变量
syms K R_p real positive

% 定义K_target
K_target = 1/(R_a - R + R_p + R_s); % Asymm.
% K_target = R^2/R_p^2 + 2*R/R_p + 4*R_s/R_p;% Symm.


% 解K_target关于R_p(K)
R_s_solution = solve(K_target - K, R_s);

% 假设R_p_solution是求解得到的R_p作为K的函数（如果可能的话）

% 定义原始表达式 y_original
y_original = wn;

% 用解得的R_p(K)替换y_original中的R_p
y_substituted = subs(y_original, R_s, R_s_solution);

% 化简
y_simplified = simplify(y_substituted);

pretty(y_simplified);