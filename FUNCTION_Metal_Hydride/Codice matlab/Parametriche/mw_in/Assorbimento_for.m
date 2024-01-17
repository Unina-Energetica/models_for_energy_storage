function dyadta = Assorbimento_4_eq(ta, ya, m_H2_in, Twa, DTml)

% Variables
m_H2_a = ya(1);     % Hydrogen mass absorbed [kg]
m_MH_a = ya(2);     % Metal hydride formed   [kg]
T_a    = ya(3);     % System's temperature   [K]

% Parameters
Ca    = 59.2;        % Kinetic constant                   [1/s]
Ea    = 21170;       % Absorption activatio nenergy       [J/mol]
R     = 8.314;       % Gas constant                       [J/mol K]
Pa    = 5E+05;       % Absorption pressure                [Pa]
m_s   = 143;         % Solid mass                         [kg]
M_H2  = 2.016/1000;  % Molar mass H2                      [kg/mol]
SC    = 3;           % Stoichiometric coefficient         [molH2/molMH]
M_MH  = 0.432;       % Molar mass MH                      [kg/mol]
Cp_H2 = 14300;       % Specific heat H2                   [J/kg K]
Cp_s  = 355;         % Specific heat solid                [J/kg K] 
Tin_a = 290;         % Inlet temperature H2               [K]
A     = 5.4;         % Area "heat interchange"            [m2]
U     = 243;         % Overall heat transfer coefficient  [W/m2 K]
DHa   = -30478;      % Absorption enthalpy                [J/mol]
DSa   = -108;        % Absorption entropy                 [J/mol K]
sl    = 0.13;        % Slope coefficient                  [-]
P0    = 1E+05;       % Reference pressure                 [Pa]

% Auxiliary equations
Peq_a = exp((DHa/(R*T_a))-DSa/R+sl*((m_MH_a/m_s)-0.5))*P0;   % Absorption equilibrium pressure [Pa]
r_a   = Ca*exp(-Ea/(R*T_a))*log(Pa/Peq_a)*(1-(m_MH_a/m_s));  % Reaction rate                   [1/s]

% ODE
d_m_H2_a_dta = m_H2_in-r_a*m_s*((M_H2*SC)/M_MH);
d_m_MH_a_dta = r_a*m_s;
d_T_a_dta    = ((m_H2_in)*Cp_H2*(Tin_a-T_a)-(A*U*DTml)-DHa*r_a*m_s*SC/M_MH)/(m_H2_a*Cp_H2+m_s*Cp_s); %Absorption

% Variables' array
dyadta = [d_m_H2_a_dta; d_m_MH_a_dta; d_T_a_dta];
end