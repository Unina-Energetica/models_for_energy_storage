function dydt = Lumped_Model_caso1des (t, y)
%Modello a parametri concentrati per il desorbimento
%% Variabili
m_H2_des = y(1);
m_MH_des = y(2);
T_des = y(3);

%% Main parameters
Vg_des = 0.0172;        % m3
Cdes = 9.6;             % 1/s
Edes = 19420;           % J/mol
R_des = 8.314;          % J/mol K
P_des = 6E+05;          % Pa
m_s_des = 143;          % kg
M_H2_des = 2.016/1000;  % kg/mol
SC_des = 3;             % -
M_MH_des = 0.432;       % kg/mol
Cp_H2_des = 14300;      % J/kg K
Cp_s_des = 355;         % J/kg K
Twd = 353;              % K
DHd = 30800;            % J/mol
DSd = 108;              % J/mol K
sl_des = 0.13;          % -
P0_des = 1E+05;         % Pa
A_des = 5.4;            % m2
U_des = 243;            % W/m2 K
m_H2_out = 2/3600;      % kg/s

%% Equazioni preliminari
Peq_des = (exp(-(DHd/(R_des*T_des))+(DSd/R_des)+sl_des*((m_MH_des/m_s_des)-0.5)))*P0_des;   % Pa
r_des = Cdes*exp(-Edes/(R_des*T_des))*((P_des-Peq_des)/Peq_des)*(m_MH_des/m_s_des);         % 1/s

%% Equazioni di governo
d_m_H2_des_dt = -m_H2_out-r_des*m_s_des*(M_H2_des*SC_des/M_MH_des);                                             % kg
d_m_MH_des_dt = r_des*m_s_des;                                                                                  % kg
d_T_des_dt = (A_des*U_des*(Twd-T_des)+DHd*r_des*m_s_des*SC_des/M_MH_des)/(m_H2_des*Cp_H2_des+m_s_des*Cp_s_des); % K

%% Variabili di stato
dydt = [d_m_H2_des_dt; d_m_MH_des_dt; d_T_des_dt];
end

