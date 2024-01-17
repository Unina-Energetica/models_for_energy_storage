clear all; clc;

%% ASSORBIMENTO 10 bar
% Parametri
Vg_a = 6.16E-05;    % Gas volume                        [m3]
Ca = 75;            % Kinetic constant                  [1/s]
Ea = 21170;         % Activation energy                 [J/mol]
R = 8.314;          % Gas constant                      [J/ mol K]
Pa = 10E+05;        % Absorption pressure               [Pa]
m_s = 500/1000;     % Solid mass                        [kg]
M_H2 = 2.016/1000;  % Molar mass                        [kg/mol]
SC = 3;             % Stoichiometric coefficient        [molH2/molMH]
M_MH = 0.421;       % Molar mass                        [kg/mol]
Cp_H2 = 14300;      % Specific heat                     [J/kg K]
Cp_s = 355;         % Specific heat                     [J/kg K] 
Tin_a = 313;        % Input temperature                 [K]
A = 3.82E-02;       % Area                              [m2]
U = 1000;           % Overall heat transfer coefficient [W/m2 K]
DHa = -28000;       % Enthalpy                          [J/mol]
DSa = -107.2;       % Entropy                           [J/mol K]
sl = 0.13;          % Slope coefficient                 [-]
P0 = 1E+05;         % Reference pressure                [Pa]
Twa = 298;          % Inlet water temperature           [K]
m_H2_in = 0.001299; % Hydrogen flow rate                [kg/s]
phi_s = 0.35;       % slope factor
phi_0 = 0.15;       % constant
phi =0.2;           % hysteresis factor

% Condizioni iniziali
    m_H2_0 = 0;
    m_MH_0 = 0;
    T0 = Twa;
    
    y0 = [m_H2_0; m_MH_0; T0];
    
    % Tempo di simulazione 
    t_in = 0;
    t_fine = 600; 
    tspan = [t_in, t_fine];
    
    ode_fun = @(t1, y1) Abs_Muthukumar(t1, y1, m_H2_in);
    
    [t1, y1] = ode45(ode_fun, tspan, y0);
    
    % Estrazione delle variabili
    m_H2_1 = y1(:, 1); % kg
    m_MH_1 = y1(:, 2); % kg
    T_1 = y1(:, 3);    % K
    
    % Equazioni addizionali
    Peq_1 = exp((DHa/(R*T_1))-DSa/R+(phi_s+phi_0)*tan(pi*((m_MH_1/m_s)-0.5)+phi/2))*P0;     % Reaction rate
    Pg_1 = abs(m_H2_1.*R.*T_1./(Vg_a.*M_H2));                           % Pa
    wt_abs_1 = (m_MH_1./m_s)*(M_H2*SC/M_MH)*100;                            % [%]
    H_M_1 = (wt_abs_1*(M_MH/M_H2))./(100-wt_abs_1);

%% ASSORBIMENTO 20 bar
Pa = 20E+05;        % Absorption pressure               [Pa]

    ode_fun = @(t2, y2) Abs_Muthukumar_20(t2, y2, m_H2_in);
    
    [t2, y2] = ode45(ode_fun, tspan, y0);
    
    % Estrazione delle variabili
    m_H2_2 = y2(:, 1); % kg
    m_MH_2 = y2(:, 2); % kg
    T_2 = y2(:, 3);    % K

    Peq_2 = exp((DHa/(R*T_2))-DSa/R+(phi_s+phi_0)*tan(pi*((m_MH_2/m_s)-0.5)+phi/2))*P0;     % Reaction rate
    Pg_2 = abs(m_H2_2.*R.*T_2./(Vg_a.*M_H2));                                               % Pa
    wt_abs_2 = (m_MH_2./m_s)*(M_H2*SC/M_MH)*100;                                            % [%]
    H_M_2 = (wt_abs_2*(M_MH/M_H2))./(100-wt_abs_2);

%% ASSORBIMENTO 30 bar
% Parametri
Pa = 30E+05;        % Absorption pressure [Pa]

    ode_fun = @(t3, y3) Abs_Muthukumar_30(t3, y3, m_H2_in);
    
    [t3, y3] = ode45(ode_fun, tspan, y0);
    
    % Estrazione delle variabili
    m_H2_3 = y3(:, 1); % kg
    m_MH_3 = y3(:, 2); % kg
    T_3 = y3(:, 3);    % K

        % Equazioni addizionali
    Peq_3 = exp((DHa/(R*T_3))-DSa/R+(phi_s+phi_0)*tan(pi*((m_MH_3/m_s)-0.5)+phi/2))*P0;     % Reaction rate
    Pg_3 = abs(m_H2_3.*R.*T_3./(Vg_a.*M_H2));                                               % Pa
    wt_abs_3 = (m_MH_3./m_s)*(M_H2*SC/M_MH)*100;                                            % [%]
    H_M_3 = (wt_abs_3*(M_MH/M_H2))./(100-wt_abs_3);
    r_3 = Ca.*exp(-Ea./(R*T_3)).*log(Pa./Peq_3).*(1-(m_MH_3./m_s));                         %Reaction rate [1/s]
  

%% Eestrapolazione dei risultati dal paper
ten_bars = readtable("10bar.xlsx");
t_10bar = ten_bars.Var1;
wt_10bar = ten_bars.Var2;

twenty_bars = readtable("20bar.xlsx");
t_20bar = twenty_bars.Var1;
wt_20bar = twenty_bars.Var2;

thirty_bars = readtable("30bar.xlsx");
t_30bar = thirty_bars.Var1;
wt_30bar = thirty_bars.Var2;

ten_Muthu = readtable("10bar_Muthukumar.xlsx");
t_10bar_M = ten_Muthu.Var1;
wt_10bar_M = ten_Muthu.Var2;

twenty_Muthu = readtable("20bar_Muthukumar.xlsx");
t_20bar_M = twenty_Muthu.Var1;
wt_20bar_M = twenty_Muthu.Var2;

thirty_Muthu = readtable("30bar_Muthukumar.xlsx");
t_30bar_M = thirty_Muthu.Var1;
wt_30bar_M = thirty_Muthu.Var2;


%% Interpolazione
wt_abs_1_interp = interp1(t1, wt_abs_1, t2);
wt_abs_3_interp = interp1(t3, wt_abs_3, t2);
wt_10bar_paper = interp1(t_10bar, wt_10bar, t2);
wt_20bar_paper = interp1(t_20bar, wt_20bar, t2);
wt_30bar_paper = interp1(t_30bar, wt_30bar, t2);
wt_10bar_Muthu = interp1(t_10bar_M, wt_10bar_M, t2);
wt_20bar_Muthu = interp1(t_20bar_M, wt_20bar_M, t2);
wt_30bar_Muthu = interp1(t_30bar_M, wt_30bar_M, t2);

%% Plot dei risultati per confronto con Talaganis e Muthukumar 
figure(1)
plot(t2, wt_abs_1_interp, 'k-', t2, wt_abs_2, 'k-', t2, wt_abs_3_interp, 'k-',t2, wt_10bar_paper, 'ro', t2, wt_20bar_paper, 'b^', t2, wt_30bar_paper, 'cs', t2, wt_10bar_Muthu, 'r--', t2, wt_20bar_Muthu, 'b--', t2, wt_30bar_Muthu, 'c--');
ylabel('Hydrogen storage capacity wt [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([0 2]);
hold on 
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
legend('Present work 10 - 20 - 30 [bar]','','','Talaganis 10 [bar]', 'Talaganis 20 [bar]','Talaganis 30 [bar]', 'Muthukumar 10 [bar]', 'Muthukumar 20 [bar]', 'Muthukumar 30 [bar]','location','northeast');
grid on
title('Effect of supply pressure on hydrogen storage capacity for MmNi_{4.6}Al_{0.4}');

%% Interpolazione
Peq_1_interp = interp1(t1, Peq_1, t2);
Peq_3_interp = interp1(t3, Peq_3, t2);
T_1_interp = interp1(t1, T_1, t2);
T_3_interp = interp1(t3, T_3, t2);
H_M_1_interp = interp1(t1, H_M_1, t2);
H_M_3_interp = interp1(t3, H_M_3, t2);

%% Plot dei rimanenti risultati per confronto con Muthukumar
figure(2)
subplot(1,1,1);
plot(t2, T_1_interp-273, 'r-', t2, T_2-273, 'b-', t2, T_3_interp-273, 'c-');
ylabel('[°C]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
ylim([0 100]);
yticks (0:5:100);
hold on
title('Average bed temperature "T"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
xlim([0 600]);
hold off
grid on
legend('10 [bar]', '20 [bar]', '30 [bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10)




% subplot(2,1,2);
figure(3)
plot(H_M_1_interp, Peq_1_interp/1E+05, 'r-', H_M_2, Peq_2/1E+05, 'b-', H_M_3_interp, Peq_3_interp/1E+05, 'c-');
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
% ylim([0 40]);
% yticks (0:5:40);
title('P-C-T characteristics', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
hold on
xlabel('Concentration (H/M)', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
% xlim([0 1.2]);
grid on
hold off

subplot(3,1,3);
plot(t3, H_M_1);
ylabel('[H/M]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
% ylim([0 40]);
% yticks (0:5:40);
title('Concentration (H/M)', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
hold on
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',10);
% xlim([0 1.2]);
grid on
hold off