%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-09-09(yyyy-mm-dd)
% ross constant-resistance T-Coil
% 1, for constant resistance
% Symmetrical Type c
%--------------------------------------------------------------------------
ConstR_TEST = 0;
syms C_b L_a L_b L_m C R s w_n zeta b R_p R_a R_s R_c R_a
AA = R_c; BB = R_a; CCC = R_a; DD = R_p; EE = R; 
ZG = R_p*(CCC*AA+EE*AA+EE*BB+EE*CCC)/(CCC*AA+CCC*BB+DD*AA+DD*BB+DD*CCC+EE*AA+EE*BB+EE*CCC);
A  = R_c/(1+s*R_c*C_b);
B  = s*L_a+R_a;
CC = s*L_a+R_a;
D  = s*L_m+R_p/(1+s*C*R_p);
E  = R;
Z_D = R_p/(1+s*C*R_p);
% ZG  = R;
CellName = {     zeta,  R,     C, R_p};
CellVal  = [sqrt(2)/2, 50, 4e-12, 500];

% Calculate constant R
Z = A*B*(CC+D+E)+A*CC*(D-E)-E^2*(A+B+CC);
FF0 = collect(Z,'s');
[nn, dd] = numden(FF0);
FF1 = collect(nn,s);
coef = coeffs(FF1, s);
Hff1 = collect(factor(coef(1)));
Hff2 = collect(factor(coef(2)));
Hff3 = collect(factor(coef(3)));
Hff4 = collect(factor(coef(4)));
y = solve(Hff1(end)==0, Hff2(end)==0, Hff3(end)==0, Hff4(end)==0, L_m, L_a, R_c, R_a);
kk = 2;
L_a0 = simplify(y.L_a(kk))% L_a = (C*R^2)/2
R_c0 = simplify(y.R_c(kk))% R_c = 4*R_p
L_m0 = simplify(y.L_m(kk))% L_m = -(R^2*(C - 4*C_b))/4
R_a0 = simplify(y.R_a(kk))% R_a = R^2/(2*R_p)
% test const R
if ConstR_TEST
    L_a_v0 = subs(L_a0, CellName, CellVal);
    R_c_v0 = subs(R_c0, CellName, CellVal);
    L_m_v0 = subs(L_m0, CellName, CellVal);
%     R_b_v0 = subs(R_b0, CellName, CellVal);
    fprintf('La=%0.3f nH\n', vpa(L_a_v0)*1e9);
    fprintf('Rc=%0.3f Ohm\n', vpa(R_c_v0));
    % fprintf('Lm=%0.3f nH\n', vpa(L_m_v0)*1e9);
%     fprintf('Rb=%0.3f Ohm\n', vpa(R_b_v0));
    return;
end

% cancelation zeros
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
y = solve(Hf2(end-1)==0, Hf3(end-1)==0,Hf5(end-1)==0,...
    L_a0==L_a,R_c0==R_c,L_m0==L_m, R_a0==R_a, R_a, L_a, R_c, L_m, w_n, C_b);
kk = 1;
L_a1  = simplify(y.L_a(kk))% L_a = (C*R^2)/2
R_c1  = simplify(y.R_c(kk))% R_c = 4*R_p
L_m1  = simplify(y.L_m(kk))% L_m = -(C*((- R^4 - 4*R^3*R_p - 4*R^2*R_p^2 + 16*R_p^4)^(1/2) + R^2 - 4*R_p^2))/4
R_a1  = simplify(y.R_a(kk))% R_a = R^2/(2*R_p)
w_n1  = simplify(y.w_n(kk))% w_n = (2^(1/2)*(2*R*R_p + (- R^4 - 4*R^3*R_p - 4*R^2*R_p^2 + 16*R_p^4)^(1/2) + R^2 + 4*R_p^2))/(2*C*R*R_p*(R + 2*R_p))
C_b1  = simplify(y.C_b(kk))% C_b = -(C*((- R^4 - 4*R^3*R_p - 4*R^2*R_p^2 + 16*R_p^4)^(1/2) - 4*R_p^2))/(4*R^2)
L_a_v = subs(L_a1, CellName, CellVal);
R_a_v = subs(R_a1, CellName, CellVal);
L_m_v = subs(L_m1, CellName, CellVal);
R_c_v = subs(R_c1, CellName, CellVal);
C_b_v = subs(C_b1, CellName, CellVal);
w_n_v = subs(w_n1, CellName, CellVal);
zeta_v = CellVal(1);
w_n_x = w_n_v*sqrt(sqrt(1+(2*zeta_v^2-1)^2)-(2*zeta_v^2-1));
fprintf('Ra=%0.3f Ohm\n', vpa(R_a_v));
fprintf('Rb=%0.3f Ohm\n', vpa(R_a_v));
fprintf('Rc=%0.3f Ohm\n', vpa(R_c_v));
fprintf('La=%0.3f nH\n', vpa(L_a_v)*1e9);
fprintf('Lb=%0.3f nH\n', vpa(L_a_v)*1e9);
fprintf('Lm=%0.3f nH\n', vpa(L_m_v)*1e9);
% fprintf('Rb=%0.3f Ohm\n', vpa(R_b_v));
fprintf('Cb=%0.3fpF(%0.4f*C)\n', vpa(C_b_v)*1e12, C_b_v/CellVal(3));
fprintf('f3dB=%0.3fGHz(%0.3f/RC)\n', vpa(w_n_x)/2/pi*1e-9, w_n_x*CellVal(2)*CellVal(3));

