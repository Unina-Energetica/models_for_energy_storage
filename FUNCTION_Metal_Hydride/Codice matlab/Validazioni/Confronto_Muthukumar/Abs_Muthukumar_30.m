function dy3dt3 = Abs_Laurencelle(t3, y3, m_H2_in)
% Variabili
m_H2_3 = y3(1);
m_MH_3 = y3(2);
T_3 = y3(3);

% Parametri
Ca = 75;            % Kinetic constant                  [1/s]
Ea = 21170;         % Activation energy                 [J/mol]
R = 8.314;          % Gas constant                      [J/ mol K]
Pa = 30E+05;        % Absorption pressure               [Pa]
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
DSa = -107.2;       % Entropy                           [J/mol K]
sl = 0.13;          % Slope coefficient                 [-]
P0 = 1E+05;         % Reference pressure                [Pa]
Twa = 298;          % Inlet water temperature           [K]
phi_s = 0.35;       % slope factor
phi_0 = 0.15;       % constant
phi =0.2;           % hysteresis factor

% Equazioni ausiliarie
Peq_3 = exp((DHa/(R*T_3))-DSa/R+(phi_s+phi_0)*tan(pi*((m_MH_3/m_s)-0.5)+phi/2))*P0;
r_3 = Ca*exp(-Ea/(R*T_3))*log(Pa/Peq_3)*(1-(m_MH_3/m_s));

% ODE
d_m_H2_3_dt3 = m_H2_in-r_3*m_s*((M_H2*SC)/M_MH);
d_m_MH_3_dt3 = r_3*m_s;
d_T_3_dt3 = ((m_H2_in)*Cp_H2*(Tin_a-T_3)+A*U*(Twa-T_3)-DHa*r_3*m_s*SC/M_MH)/(m_H2_3*Cp_H2+m_s*Cp_s);

% Variabili di stato in vettore colonna
dy3dt3 = [d_m_H2_3_dt3; d_m_MH_3_dt3; d_T_3_dt3];
end