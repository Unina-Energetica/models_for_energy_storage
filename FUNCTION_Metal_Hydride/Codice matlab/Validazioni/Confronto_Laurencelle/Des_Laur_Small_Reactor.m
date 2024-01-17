function dyddtd = Des_Laurencelle(td, yd, m_H2_out)
% Variabili
m_H2_d = yd(1);
m_MH_d = yd(2);
T_d = yd(3);

% Parametri
Cd = 9.6;          % 1/s
Ed = 16420;         % J/mol
R = 8.314;          % J/ mol K
Pd = 0.068E+05;         % Pa
m_s = 1/1000;          % kg
M_H2 = 2.016/1000;  % kg/mol
SC = 2.76;             % molH2/molMH
M_MH = 0.432;       % kg/mol
Cp_H2 = 14300;      % J/kg K
Cp_s = 355;         % J/kg K 
A = 5E-04;            % m2
U = 80;            % W/m2 K
DHd = 30800;       % J/mol
DSd = 108;         % J/mol K
sl = 0.13;          % -
P0 = 1E+05;         % Pa
Twd = 296;

% Equazioni ausiliarie
Peq_d = (exp(-(DHd/(R*T_d))+(DSd/R)+sl*((m_MH_d/m_s)-0.5)))*P0;
r_d = Cd*exp(-Ed/(R*T_d))*((Pd-Peq_d)/Peq_d)*(m_MH_d/m_s);

% ODE
d_m_H2_d_dtd = -m_H2_out-r_d*m_s*((M_H2*SC)/M_MH);
d_m_MH_d_dtd = r_d*m_s;
d_T_d_dtd = (A*U*(Twd-T_d)+DHd*r_d*m_s*SC/M_MH)/(m_H2_d*Cp_H2+m_s*Cp_s);

% Variabili di stato in vettore colonna
dyddtd = [d_m_H2_d_dtd; d_m_MH_d_dtd; d_T_d_dtd];
end