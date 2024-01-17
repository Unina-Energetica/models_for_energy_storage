function dyadta = Abs_Laurencelle(ta, ya, m_H2_in)
% Variabili
m_H2_a = ya(1);
m_MH_a = ya(2);
T_a = ya(3);

% Parametri
Ca = 59.2;          % 1/s
Ea = 21170;         % J/mol
R = 8.314;          % J/ mol K
Pa = 6E+05;         % Pa
m_s = 1/1000;       % kg
M_H2 = 2.016/1000;  % kg/mol
SC = 2.76;          % molH2/molMH
M_MH = 0.432;       % kg/mol
Cp_H2 = 14300;      % J/kg K
Cp_s = 355;         % J/kg K 
Tin_a = 290;        % K
A = 5E-04;          % m2
U = 80;             % W/m2 K
DHa = -30478;       % J/mol
DSa = -108;         % J/mol K
sl = 0.13;          % -
P0 = 1E+05;         % Pa
Twa = 296;          % K

% Equazioni ausiliarie
Peq_a = exp((DHa/(R*T_a))-DSa/R+sl*((m_MH_a/m_s)-0.5))*P0;
r_a = Ca*exp(-Ea/(R*T_a))*log(Pa/Peq_a)*(1-(m_MH_a/m_s));

% ODE
d_m_H2_a_dta = m_H2_in-r_a*m_s*((M_H2*SC)/M_MH);
d_m_MH_a_dta = r_a*m_s;
d_T_a_dta = ((m_H2_in)*Cp_H2*(Tin_a-T_a)+A*U*(Twa-T_a)-DHa*r_a*m_s*SC/M_MH)/(m_H2_a*Cp_H2+m_s*Cp_s);

% Variabili di stato in vettore colonna
dyadta = [d_m_H2_a_dta; d_m_MH_a_dta; d_T_a_dta];
end