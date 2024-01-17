function dyddtd = Desorbimento(td, yd, m_H2_out, Twd)
% Variables
m_H2_d = yd(1);     % Hydrogen mass desorbed [kg]
m_MH_d = yd(2);     % Metal hydride left     [kg]
T_d = yd(3);        % System's temperature   [K]

% Parameters
Cd = 9.6;           % Kinetic constant                   [1/s]
Ed = 16420;         % Desorption activatio nenergy       [J/mol]
R = 8.314;          % Gas constant                       [J/ mol K]
Pd = 6E+05;         % Desorption pressure                [Pa]
m_s = 143;          % Solid mass                         [kg]
M_H2 = 2.016/1000;  % Molar mass H2                      [kg/mol]
SC = 3;             % Stoichiometric coefficient         [molH2/molMH]
M_MH = 0.432;       % Molar mass MH                      [kg/mol]
Cp_H2 = 14300;      % Specific heat H2                   [J/kg K]
Cp_s = 355;         % Specific heat solid                [J/kg K] 
A = 5.4;            % Area "heat interchange"            [m2]
U = 243;            % Overall heat transfer coefficient  [W/m2 K]
DHd = 30800;        % Absorption enthalpy                [J/mol]
DSd = 108;          % Absorption entropy                 [J/mol K]
sl = 0.13;          % Slope coefficient                  [-]
P0 = 1E+05;         % Reference pressure                 [Pa]

% Auxiliary equations
Peq_d = (exp(-(DHd/(R*T_d))+(DSd/R)+sl*((m_MH_d/m_s)-0.5)))*P0; % Desorption equilibrium pressure [Pa]
r_d = Cd*exp(-Ed/(R*T_d))*((Pd-Peq_d)/Peq_d)*(m_MH_d/m_s);      % Reaction rate

% ODE
d_m_H2_d_dtd = -m_H2_out-r_d*m_s*((M_H2*SC)/M_MH);
d_m_MH_d_dtd = r_d*m_s;
d_T_d_dtd = (A*U*(Twd-T_d)+DHd*r_d*m_s*SC/M_MH)/(m_H2_d*Cp_H2+m_s*Cp_s);

% Variables' array
dyddtd = [d_m_H2_d_dtd; d_m_MH_d_dtd; d_T_d_dtd];
end