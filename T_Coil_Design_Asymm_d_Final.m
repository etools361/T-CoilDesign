%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-09-10(yyyy-mm-dd)
% ross T-Coil
% 1, for constant resistance
% 2, for zeros cancellation
% Asymmetrical Type d
%--------------------------------------------------------------------------
ConstR_TEST = 0;
syms C_b L_a L_b L_m C R s w_n zeta b R_p R_b R_s R_a
A  = 1/(s*C_b);
B  = s*L_a+R_a;
CC = s*L_b+R_b;
D  = s*L_m+R_p/(1+s*C*R_p)+R_s;
E  = R;
Z_D = R_p/(1+s*C*R_p);
ZG  = R_p/(1+(R_s+R_p)/(R_b+R));
% define test value
CellName = {     zeta,  R,     C, R_p, R_s,  R_a};
CellVal  = [sqrt(2)/2, 50, 4e-12, 500,  10, 1.907];

% constant resistance
Z = A*B*(CC+D+E)+A*CC*(D-E)-E^2*(A+B+CC);
FF0 = collect(Z,'s');
[nn, dd] = numden(FF0);
FF1 = collect(nn,s);
coef = coeffs(FF1, s);
y = solve(coef(1)==0, coef(2)==0, coef(3)==0,coef(4)==0, L_b, L_m, L_a, R_b);
kk = 1;% select index
L_a0 = simplify(y.L_a(kk))
L_b0 = simplify(y.L_b(kk))
L_m0 = simplify(y.L_m(kk))
R_b0 = simplify(y.R_b(kk))
% test const R
if ConstR_TEST
    L_a_v0 = subs(L_a0, CellName, CellVal);
    L_b_v0 = subs(L_b0, CellName, CellVal);
    L_m_v0 = subs(L_m0, CellName, CellVal);
    R_b_v0 = subs(R_b0, CellName, CellVal);
    fprintf('La=%0.3f nH\n', vpa(L_a_v0)*1e9);
    fprintf('Lb=%0.3f nH\n', vpa(L_b_v0)*1e9);
    % fprintf('Lm=%0.3f nH\n', vpa(L_m_v0)*1e9);
    fprintf('Rb=%0.3f Ohm\n', vpa(R_b_v0));
    return;
end

% zeros cancellation
H = Z_D*(CC*A+E*A+E*B+E*CC)/(CC*A+CC*B+D*A+D*B+D*CC+E*A+E*B+E*CC);
H0 = collect(H,'s');
[Hn, Hd] = numden(H0);
Hz = Hd*w_n^2*ZG - Hn*(s^2+2*zeta*w_n*s+w_n^2);
Hz0 = collect(Hz,s);
Hcoef = coeffs(Hz0, s);
Hf1 = collect(factor(Hcoef(1)));
Hf2 = collect(factor(Hcoef(2)));
Hf3 = collect(factor(Hcoef(3)));
Hf4 = collect(factor(Hcoef(4)));
Hf5 = collect(factor(Hcoef(5)));
y = solve(Hf2(end-1)==0, Hf3(end-1)==0,Hf4(end-1)==0,Hf5(end-1)==0,...
    L_a0==L_a,L_b0==L_b,L_m0==L_m, R_b0 == R_b, R_b, L_a, L_b, L_m, w_n, C_b);
kk = 1;% select index
L_a1  = simplify(y.L_a(kk))
L_b1  = simplify(y.L_b(kk))
L_m1  = simplify(y.L_m(kk))
R_b1  = simplify(y.R_b(kk))
w_n1  = simplify(y.w_n(kk))
C_b1  = simplify(y.C_b(kk))

L_a_v = subs(L_a1, CellName, CellVal);
L_b_v = subs(L_b1, CellName, CellVal);
L_m_v = subs(L_m1, CellName, CellVal);
R_b_v = subs(R_b1, CellName, CellVal);
C_b_v = subs(C_b1, CellName, CellVal);
w_n_v = subs(w_n1, CellName, CellVal);
zeta_v = CellVal(1);
w_n_x = w_n_v*sqrt(sqrt(1+(2*zeta_v^2-1)^2)-(2*zeta_v^2-1));
fprintf('La=%0.3f nH\n', vpa(L_a_v)*1e9);
fprintf('Lb=%0.3f nH\n', vpa(L_b_v)*1e9);
fprintf('Lm=%0.3f nH\n', vpa(L_m_v)*1e9);
fprintf('Rb=%0.3f Ohm\n', vpa(R_b_v));
fprintf('Cb=%0.3fpF(%0.4f*C)\n', vpa(C_b_v)*1e12, C_b_v/CellVal(3));
fprintf('f3dB=%0.3fGHz(%0.3f/RC)\n', vpa(w_n_x)/2/pi*1e-9, w_n_x*CellVal(2)*CellVal(3));
Gain0=simplify(subs(ZG, {R_b}, {R_b1}));
Gain = subs(Gain0/R, CellName, CellVal);
fprintf('Gain=%0.3f\n', vpa(Gain));

% 验证
fprintf('---Asymm.验证公式---\n');
Ra = CellVal(6); Rp = CellVal(4); Rs = CellVal(5); R = CellVal(2); C = CellVal(3); zeta = CellVal(1);
T2 = 1 / (Ra + Rp + Rs - R);
T3 = T2 * (R - Ra) + 1;
T4 = (T2 * Rp) / (T3^2);
Rb = T2 * (R - Ra)^2 - Ra;
La = Rp * C * (R - Ra) * (sqrt(Rp * T2) - 1);
Lb = Rp * C * (R - Ra) * (T3 - sqrt(Rp * T2));
Cb = C * (Rp * T3) / (R * (T3-1)) * (2 * zeta^2 + sqrt(T4) - 1 - 2 * zeta * sqrt(zeta^2 + sqrt(T4) - 1));
Lm = R^2 * Cb - (La * Lb) / (La + Lb);
wn = -( zeta + sqrt(zeta^2 + sqrt(T4) - 1)) / (C * Rp * (sqrt(T4) - 1));
w_n_x = wn*sqrt(sqrt(1+(2*zeta^2-1)^2)-(2*zeta^2-1));

fprintf('La=%0.3f nH\n', La*1e9);
fprintf('Lb=%0.3f nH\n', Lb*1e9);
fprintf('Lm=%0.3f nH\n', Lm*1e9);
fprintf('Rb=%0.3f Ohm\n', Rb);
fprintf('Cb=%0.3fpF(%0.4f*C)\n', Cb*1e12, Cb/CellVal(3));
fprintf('f3dB=%0.3fGHz(%0.3f/RC)\n', vpa(w_n_x)/2/pi*1e-9, w_n_x*CellVal(2)*CellVal(3));
Gain0 = ((1 - R_a/R))/(1 + R_s/R_p);
Gain = subs(Gain0, CellName, CellVal);
fprintf('Gain=%0.3f\n', vpa(Gain));
