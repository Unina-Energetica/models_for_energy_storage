function dyabsdtabs = Lumped_Model(tabs, yabs, m_H2_in)
%Modello nuovo a parametri concentrati
%Variabili
m_H2_abs = yabs(1);
m_MH_abs = yabs(2);
T_abs = yabs(3);

%Main parameters
Cabs = 59.2;          % 1/s
Eabs = 21170;         % J/mol
R_abs = 8.314;        % J/ mol K
P_abs = 5E+05;        % Pa
m_s_abs = 143;        % kg
M_H2_a = 2.016/1000;  % kg/mol
SC_a = 3;             % molH2/molMH
M_MH_a = 0.432;       % kg/mol
Cp_H2_abs = 14300;    % J/kg K
Cp_s_abs = 355;       % J/kg K 
Tin_abs = 290;        % K
A_abs = 5.4;          % m2
U_abs = 243;          % W/m2 K
Twa = 298;            % K
DHa = -30478;         % J/mol
DSa = -108;           % J/mol K
sl_abs = 0.13;        % -
P0_abs = 1E+05;       % Pa

Peq_abs = exp((DHa/(R_abs*T_abs))-DSa/R_abs+sl_abs*((m_MH_abs/m_s_abs)-0.5))*P0_abs;
r_abs = Cabs*exp(-Eabs/(R_abs*T_abs))*log(P_abs/Peq_abs)*(1-(m_MH_abs/m_s_abs));

%ODE
d_m_H2_abs_dtabs = m_H2_in-r_abs*m_s_abs*((M_H2_a*SC_a)/M_MH_a);
d_m_MH_abs_dtabs = r_abs*m_s_abs;
d_T_abs_dtabs = ((m_H2_in)*Cp_H2_abs*(Tin_abs-T_abs)+A_abs*U_abs*(Twa-T_abs)-DHa*r_abs*m_s_abs*SC_a/M_MH_a)/(m_H2_abs*Cp_H2_abs+m_s_abs*Cp_s_abs);

%Variabili di stato in vettore colonna
dyabsdtabs = [d_m_H2_abs_dtabs; d_m_MH_abs_dtabs; d_T_abs_dtabs];
end


