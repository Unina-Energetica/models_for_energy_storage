function dy2dt2 = Abs_Laurencelle(t2, y2, m_H2_in)
% Variabili
m_H2_2 = y2(1);
m_MH_2 = y2(2);
T_2 = y2(3);

% Parametri
Ca = 75;            % Kinetic constant                  [1/s]
Ea = 21170;         % Activation energy                 [J/mol]
R = 8.314;          % Gas constant                      [J/ mol K]
Pa = 20E+05;        % Absorption pressure               [Pa]
m_s = 500/1000;     % Solid mass                        [kg]
M_H2 = 2.016/1000;  % Molar mass                        [kg/mol]
SC = 3;             % Stoichiometric coefficient        [molH2/molMH]
M_MH = 0.421;       % Molar mass                        [kg/mol]
Cp_H2 = 14300;      % Specific heat                     [J/kg K]
Cp_s = 355;         % Specific heat                     [J/kg K] 
Tin_a = 313;        % Input temperature                 [K]
A = 3.82E-02;       % Area                              [m2]
U = 1000;            % Overall heat transfer coefficient [W/m2 K]
DHa = -28000;       % Enthalpy                          [J/mol]
DSa = -107.2;       % Entroyp                           [J/mol K]
sl = 0.13;          % Slope coefficient                 [-]
P0 = 1E+05;         % Reference pressure                [Pa]
Twa = 298;          % Inlet water temperature           [K]
phi_s = 0.35;       % slope factor
phi_0 = 0.15;       % constant
phi =0.2;           % hysteresis factor

% Equazioni ausiliarie
Peq_2 = exp((DHa/(R*T_2))-DSa/R+(phi_s+phi_0)*tan(pi*((m_MH_2/m_s)-0.5)+phi/2))*P0;
r_2 = Ca*exp(-Ea/(R*T_2))*log(Pa/Peq_2)*(1-(m_MH_2/m_s));

% ODE
d_m_H2_2_dt2 = m_H2_in-r_2*m_s*((M_H2*SC)/M_MH);
d_m_MH_2_dt2 = r_2*m_s;
d_T_2_dt2 = ((m_H2_in)*Cp_H2*(Tin_a-T_2)+A*U*(Twa-T_2)-DHa*r_2*m_s*SC/M_MH)/(m_H2_2*Cp_H2+m_s*Cp_s);

% Variabili di stato in vettore colonna
dy2dt2 = [d_m_H2_2_dt2; d_m_MH_2_dt2; d_T_2_dt2];
end