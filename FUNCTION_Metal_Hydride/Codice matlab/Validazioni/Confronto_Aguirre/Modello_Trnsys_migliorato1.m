clear all; clc;
tic
%% Parameters absorption and desorption
Vg    = 0.0172;      % Gas Volume                         [m3]
R     = 8.314;       % Gas constant                       [J/mol K]
M_H2  = 2.016/1000;  % Molar mass H2                      [kg/mol]
M_MH  = 0.432;       % Molar mass   MH                    [kg/mol]
Cp_H2 = 14300;       % Specific heat H2                   [J/kg K]
Cp_s  = 355;         % Specific heat solid                [J/kg K]
A     = 5.4;         % Area "heat interchange"            [m2]
U     = 243;         % Overall heat transfer coefficient  [W/m2 K]
sl    = 0.13;        % Slope coefficient                  [-]
P0    = 1E+05;       % PReference pressure                [Pa]
m_s   = 143;         % MSolid mass                        [kg]
SC    = 3;           % Stoichiometric coefficient

%% Input parameters
m_H2_in  = 2/3600;  % Input hydrogen flowrate            [kg/s]
m_H2_out = 2/3600;  % Output hydrogen flowrate           [kg/s]
Twa      = 298;     % Cooling water temperature          [K]
Twd      = 353;     % Heating water temperature          [K]

%% ABSORPTION

if (m_H2_in > 0) %&& (m_H2_out == 0) 

    %Absorption parameters
    Ca    = 59.2;   % Kinetic constant                   [1/s]
    Ea    = 21170;  % Absorption activation energy       [J/mol]
    Pa    = 5E+05;  % Absorption pressure                [Pa]
    Tin_a = 290;    % H2 input temperature               [K]
    DHa   = -30478; % Absorption enthalpy                [J/mol]
    DSa   = -108;   % Absorption entropy                 [J/mol K]

    %Initial conditions
    m_H20_a = 0;    % Absorbed hydrogen mass             [kg]
    m_MH0_a = 0;    % Metal hydride mass                 [kg]
    T0_a    = Twa;  % System's temperature               [K]

    y0_a    = [m_H20_a; m_MH0_a; T0_a]; % Initial conditions array

    %Simulation time
    t_start_a = 0;
    t_end_a   = 1800;

    t_span_a  = [t_start_a, t_end_a];

    %Function handle 
    ode_fun_a = @(ta, ya) Assorbimento(ta, ya, m_H2_in, Twa);

    %Differential Equations solver
    [ta, ya] = ode23s(ode_fun_a, t_span_a, y0_a);

    %Variables' extraction
    m_H2_a = ya(:,1);   % [kg]
    m_MH_a = ya(:,2);   % [kg]
    T_a    = ya(:,3);   % [K]

    %Other equations
    Peq_a = exp((DHa./(R*T_a))-DSa./R+sl.*((m_MH_a./m_s)-0.5))*P0;    % Absorption equilibrium pressure          [Pa]
    r_a   = Ca.*exp(-Ea./(R*T_a)).*log(Pa./Peq_a).*(1-(m_MH_a./m_s)); % Reaction rate                            [1/s]
    Pg_a  = (m_H2_a*R.*T_a./(Vg.*(M_H2/1000)));                       % Gas pressure                             [Pa]
    wt_a  = (m_MH_a./m_s)*(M_H2*SC/M_MH)*100;                         % Hydrogenation capacity in weight percent [%]
    H_M_a = (wt_a.*(M_MH/M_H2))./(100.-wt_a);                         % Concentration ratio                      [-]
    % Peq_poly = (96.4158.*(H_M_a)-484.7722.*(H_M_a.^2)+1117.237.*(H_M_a.^3) ...
    %     -1193.2852.*(H_M_a.^4)+479.9432.*(H_M_a.^5)).*exp(-3323.884.*((1./T_a)-(1/333)))*10E+05;
end

%% DESORPTION

if m_H2_out > 0 %elseif (m_H2_in == 0) && (m_H2_out > 0)

    %Desorption parameters 
    Cd  = 9.6;      % Kinetic constant              [1/s]
    Ed  = 16420;    % Desorption activation energy  [J/mol]
    Pd  = 6E+05;    % Desorption pressure           [Pa]
    DHd = 30800;    % Desorption enthalpy           [J/mol]
    DSd = 108;      % Desorption entropy            [J/mol K]

    %Initial conditions
    m_H20_d = 0;            % Desorbed hydrogen mass [kg]
    m_MH0_d = m_MH_a(end);  % Metal hydride left     [kg]
    T0_d    = 352.61;       % System temperature     [K]

    y0_d    = [m_H20_d; m_MH0_d; T0_d]; % Initial conditions array

    %Simulation time
    t_start_d = t_end_a+1;
    t_end_d   = 3600;

    t_span_d  = [t_start_d, t_end_d];
    

    %Function handle 
    ode_fun_d = @(td, yd) Desorbimento(td, yd, m_H2_out, Twd);

    %Differential Equations solver
    [td, yd]  = ode23s(ode_fun_d,t_span_d,y0_d);

    %Variables' extraction
    m_H2_d = yd(:,1);   % [kg]
    m_MH_d = yd(:,2);   % [kg]
    T_d    = yd(:,3);   % [K]

    %Other equations
    Peq_d = exp(-(DHd./(R.*T_d))+DSd./R+sl.*((m_MH_d./m_s)-0.5))*P0;   % Desorption equilibrium pressure             [Pa]
    r_d   = Cd.*exp(-Ed./(R.*T_d)).*((Pd-Peq_d)./Peq_d).*(m_MH_d./m_s);% Reaction rate                               [1/s]
    Pg_d  = (m_H2_d*R.*T_d./(Vg.*(M_H2/1000)));                        % Hydrogenation capacity in weight percent    [%]
    wt_d  = (m_MH_d./m_s)*(M_H2*SC/M_MH)*100;                          % Concentration ratio                         [-]
    H_M_d = (wt_d.*(M_MH/M_H2))./(100.-wt_d);                          % Hydrogen to metal ratio                     [-]
end

%% Antonio's Results plot
figure(1)
subplot(4, 1, 1);
yyaxis left
plot(ta, m_MH_a);
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
title('Metal Hydride mass "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
hold on
yyaxis right
plot(td, m_MH_d);
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
grid on

subplot(4, 1, 2);
yyaxis left
plot(ta, T_a);
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
hold on
yyaxis right
plot(td, T_d);
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('System Temperature "T"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(4, 1, 3);
yyaxis left
plot(ta, Peq_a/1E+05);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 15]);
hold on
yyaxis right
plot(td, Peq_d/1E+05);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 15]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Equilibrium pressure "P_{eq}"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(4, 1, 4);
yyaxis left
plot(ta, r_a);
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-0.0030 0.0030]);
hold on
yyaxis right
plot(td, r_d);
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-0.0030 0.0030]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Reaction rate "r"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

sgtitle('Sommella Results', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);

figure(2)
subplot(2, 2, 1);
yyaxis left
plot(ta, wt_a);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 2]);
hold on
yyaxis right
plot(td, wt_d);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 2]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Hydrogenation capacity in weight percent "wt"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(2, 2, 2);
yyaxis left
plot(ta, m_H2_a);
ylabel('[kg/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 5]);
hold on
yyaxis right
plot(td, m_H2_d);
ylabel('[kg/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 5]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Hydrogen mass absorbed "m_{H_2}"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on

subplot(2, 2, 3);
yyaxis left
plot(ta, Pg_a/1E+05);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%ylim([-10 0]);
hold on
yyaxis right
plot(td, Pg_d/1E+05);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%ylim([-5 5]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Gas pressure "P_g"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(2,2,4);
yyaxis left
plot(ta, H_M_a);
ylabel('[H/M]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 5]);
hold on 
yyaxis right
plot(td, H_M_d);
ylabel('[H/M]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 5]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Concentration ratio "H/M"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on
 
sgtitle('Other Results', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);

figure
plot(ta, T_a);
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
title('Old system Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
hold on
%% Importing the results from Talaganis, Meyer and Aguirre paper 
% m_MH [kg] absorption
m_MH_abs_paper = xlsread('m_MH_abs_copia.xlsx','Foglio1','B1:B129');
time_m_MH_abs_paper = xlsread('m_MH_abs_copia.xlsx','Foglio1','A1:A129');

% m_MH [kg] desorption
m_MH_des_paper = xlsread('m_MH_desr.xlsx','Foglio1','B1:B129');
time_m_MH_des_paper = xlsread('m_MH_desr.xlsx','Foglio1','A1:A129');


% T [K] absorption
T_abs_paper = xlsread('T_abs.xlsx','Foglio1','B1:B168');
time_T_abs_paper = xlsread('T_abs.xlsx','Foglio1','A1:A168');

% T [K] desorption
T_des_paper = xlsread('T_des.xlsx','Foglio1','B1:B181');
time_T_des_paper = xlsread('T_des.xlsx','Foglio1','A1:A181');


% Peq [bar] absorption
Peq_abs_paper = xlsread('Peq_abs.xlsx','Foglio1','B1:B178');
time_Peq_abs_paper = xlsread('Peq_abs.xlsx','Foglio1','A1:A178');

% Peq [bar] desorption
Peq_des_paper = xlsread('Peq_des.xlsx','Foglio1','B1:B203');
time_Peq_des_paper = xlsread('Peq_des.xlsx','Foglio1','A1:A203');


% reaction rate absorption
r_abs_paper = xlsread('reaction_rate_abs.xlsx','Foglio1','B1:B116');
time_r_abs_paper = xlsread('reaction_rate_abs.xlsx','Foglio1','A1:A116');

% reaction rate desorption
r_des_paper = xlsread('reaction_rate_des.xlsx','Foglio1','B1:B118');
time_r_des_paper = xlsread('reaction_rate_des.xlsx','Foglio1','A1:A118');


%% Talaganis' Results plot
figure(3)
subplot(4, 1, 1);
yyaxis left
plot(time_m_MH_abs_paper, m_MH_abs_paper);
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
hold on
yyaxis right
plot(time_m_MH_des_paper, m_MH_des_paper);
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Metal Hydride mass "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(4, 1, 2);
yyaxis left
plot(time_T_abs_paper, T_abs_paper);
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
hold on
yyaxis right
plot(time_T_des_paper, T_des_paper);
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('System Temperature "T"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

subplot(4, 1, 3);
yyaxis left
plot(time_Peq_abs_paper, Peq_abs_paper);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 15]);
hold on
yyaxis right
plot(time_Peq_des_paper, Peq_des_paper);
ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 15]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
title('Equilibrium pressure "P_{eq}"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
grid on

subplot(4, 1, 4);
yyaxis left
plot(time_r_abs_paper, r_abs_paper);
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-0.0030 0.0030]);
hold on
yyaxis right
plot(time_r_des_paper, r_des_paper);
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-0.0030 0.0030]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Reaction rate "r"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on

sgtitle('Talaganis Results','FontName','Times New Roman','FontWeight','bold', 'FontSize',18);

%% Results interpolation

% m_MH 
m_MH_a_interp = interp1(ta, m_MH_a, time_m_MH_abs_paper);
m_MH_d_interp = interp1(td, m_MH_d, time_m_MH_des_paper);

% T 
T_a_interp = interp1(ta, T_a, time_T_abs_paper);
T_d_interp = interp1(td, T_d, time_T_des_paper);

% Peq 
Peq_a_interp = interp1(ta, Peq_a, time_Peq_abs_paper);
Peq_d_interp = interp1(td, Peq_d, time_Peq_des_paper);

% r 
r_a_interp = interp1(ta, r_a, time_r_abs_paper);
r_d_interp = interp1(td, r_d, time_r_des_paper);

%% Results Comparison
figure(4)
subplot(2,2,1);
yyaxis left
plot(time_m_MH_abs_paper, m_MH_abs_paper, 'b-', time_m_MH_des_paper, m_MH_des_paper, 'b--');
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
yyaxis right
plot(time_m_MH_abs_paper, m_MH_a_interp, 'r-', time_m_MH_des_paper, m_MH_d_interp, 'r--');
ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-30 180]);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Metal Hydride mass "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on
legend('m_{MH,}_a_{,Talaganis}', 'm_{MH,}_d_{,Talaganis}', 'm_{MH,}_a_{,Antonio}', 'm_{MH,}_d_{,Antonio}');

subplot(2,2,2);
yyaxis left
plot(time_T_abs_paper, T_abs_paper, 'b-', time_T_des_paper, T_des_paper, 'b--');
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
yyaxis right
plot(time_T_abs_paper, T_a_interp, 'r-', time_T_des_paper, T_d_interp, 'r--');
ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([290 360]);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('System Temperature "T"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on
legend('T_a_{,Talaganis}', 'T_d_{,Talaganis}', 'T_a_{,Antonio}', 'T_d_{,Antonio}','location','northwest');

subplot(2,2,3);
yyaxis left
plot(time_Peq_abs_paper, Peq_abs_paper, 'b-', time_Peq_des_paper, Peq_des_paper, 'b--');
ylabel('[bar]');
ylim([0 15]);
yyaxis right
plot(time_Peq_abs_paper, Peq_a_interp/1E+05, 'r-', time_Peq_des_paper, Peq_d_interp/1E+05, 'r--');
ylabel('[bar]');
ylim([0 15]);
xlabel('Time [s]');
xlim([0 3600]);
title('Equilibrium pressure "P_{eq}"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on
legend('P_{eq,}_a_{,Talaganis}', 'P_{eq,}_d_{,Talaganis}', 'P_{eq,}_a_{,Antonio}', 'P_{eq,}_d_{,Antonio}','location','northwest');

subplot(2,2,4);
yyaxis left
plot(time_r_abs_paper, r_abs_paper, 'b-', time_r_des_paper, r_des_paper, 'b--');
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-0.0030 0.0030]);
yyaxis right
plot(time_r_abs_paper, r_a_interp, 'r-', time_r_des_paper, r_d_interp, 'r--');
ylim([-0.0030 0.0030]);
ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('Reaction rate "r"','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
grid on
legend('r_a_{,Talaganis}', 'r_d_{,Talaganis}', 'r_a_{,Antonio}', 'r_d_{,Antonio}');

sgtitle('Results Comparison','FontName','Times New Roman','FontWeight','bold', 'FontSize',18);

%% Percentage error calculation
% m_MH
err_m_MH_abs = abs(((m_MH_a_interp - m_MH_abs_paper) ./ m_MH_a_interp) * 100);
err_m_MH_des = abs(((m_MH_d_interp - m_MH_des_paper) ./ m_MH_d_interp) * 100);
nan_m_MH_des = isnan(err_m_MH_des);
err_m_MH_des(nan_m_MH_des) = 0;
mean_err_m_MH_abs = mean(err_m_MH_abs);
mean_err_m_MH_des = mean(err_m_MH_des);

% T
err_T_abs = abs(((T_a_interp - T_abs_paper) ./ T_a_interp) * 100);
err_T_des = abs(((T_d_interp - T_des_paper) ./ T_d_interp) * 100);
nan_t_des = isnan(err_T_des);
err_T_des(nan_t_des) = 0;
mean_err_T_abs = mean(err_T_abs);
mean_err_T_des = mean(err_T_des);

% Peq
err_Peq_abs = abs((((Peq_a_interp/1E+05) - Peq_abs_paper) ./ (Peq_a_interp/1E+05))* 100);
err_Peq_des = abs((((Peq_d_interp/1E+05) - Peq_des_paper) ./ (Peq_d_interp/1E+05)) * 100);
nan_Peq_des = isnan(err_Peq_des);
err_Peq_des(nan_Peq_des) = 0;
mean_err_Peq_abs = mean(err_Peq_abs);
mean_err_Peq_des = mean(err_Peq_des);

% r
err_r_abs = abs(((r_a_interp - r_abs_paper) ./ r_a_interp) * 100);
err_r_des = abs(((r_d_interp - r_des_paper) ./ r_d_interp) * 100);
nan_r_des = isnan(err_r_des);
err_r_des(nan_r_des) = 0;
mean_err_r_abs = mean(err_r_abs);
mean_err_r_des = mean(err_r_des);

%% Percentage error plots
figure(5)
subplot(2, 2, 1);
yyaxis left
plot(time_m_MH_abs_paper, err_m_MH_abs);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold on
yyaxis right
plot(time_m_MH_des_paper, err_m_MH_des);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('error "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on

subplot(2, 2, 2);
yyaxis left
plot(time_T_abs_paper, err_T_abs);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold on
yyaxis right
plot(time_T_des_paper, err_T_des);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('error "T"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on

subplot(2, 2, 3);
yyaxis left
plot(time_Peq_abs_paper, err_Peq_abs);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold on
yyaxis right
plot(time_Peq_des_paper, err_Peq_des);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('error "P_{eq}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on

subplot(2, 2, 4);
yyaxis left
plot(time_r_abs_paper, err_r_abs);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold on
yyaxis right
plot(time_r_des_paper, err_r_des);
ylabel('[%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([-5 100]);
hold off
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
xlim([0 3600]);
title('error "r"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13)
grid on

sgtitle('Percentage errors', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);

%% Normal distribution calculation
% Error m_MH absorption
[media1, mediana1, gaussiana1, deviazione_standard1] = Statistica(err_m_MH_abs);
x_m_MH_abs = linspace(media1 - 4*deviazione_standard1, media1 + 4*deviazione_standard1, 100);
y_m_MH_abs = normpdf(x_m_MH_abs, media1, deviazione_standard1);
% Error m_MH desorption
[media2, mediana2, gaussiana2, deviazione_standard2] = Statistica(err_m_MH_des);
x_m_MH_des = linspace(media2 - 4*deviazione_standard2, media2 + 4*deviazione_standard2, 100);
y_m_MH_des = normpdf(x_m_MH_des, media2, deviazione_standard2);


% Error T absorption
[media3, mediana3, gaussiana3, deviazione_standard3] = Statistica(err_T_abs);
x_T_abs = linspace(media3 - 4*deviazione_standard3, media3 + 4*deviazione_standard3, 100);
y_T_abs = normpdf(x_T_abs, media3, deviazione_standard3);
% Error T desorption
[media4, mediana4, gaussiana4, deviazione_standard4] = Statistica(err_T_des);
x_T_des = linspace(media4 - 4*deviazione_standard4, media4 + 4*deviazione_standard4, 100);
y_T_des = normpdf(x_T_des, media4, deviazione_standard4);


% Error Peq absorption
[media5, mediana5, gaussiana5, deviazione_standard5] = Statistica(err_Peq_abs);
x_Peq_abs = linspace(media5 - 4*deviazione_standard5, media5 + 4*deviazione_standard5, 100);
y_Peq_abs = normpdf(x_Peq_abs, media5, deviazione_standard5);
% Error Peq desorption
[media6, mediana6, gaussiana6, deviazione_standard6] = Statistica(err_Peq_des);
x_Peq_des = linspace(media6 - 4*deviazione_standard6, media6 + 4*deviazione_standard6, 100);
y_Peq_des = normpdf(x_Peq_des, media6, deviazione_standard6);


% Error r absorption
[media7, mediana7, gaussiana7, deviazione_standard7] = Statistica(err_r_abs);
x_r_abs = linspace(media7 - 4*deviazione_standard7, media7 + 4*deviazione_standard7, 100);
y_r_abs = normpdf(x_r_abs, media7, deviazione_standard7);
% Error r desorption
[media8, mediana8, gaussiana8, deviazione_standard8] = Statistica(err_r_des);
x_r_des = linspace(media8 - 4*deviazione_standard8, media8 + 4*deviazione_standard8, 100);
y_r_des = normpdf(x_r_des, media8, deviazione_standard8);


%% Normal distribution plots
figure(6)
subplot(4,2,1);
title('m_{MH,}_{abs}')
plot(x_m_MH_abs, y_m_MH_abs);
title('Normal distribution of m_{MH,}_{abs}');
grid on;

subplot(4,2,2);
title('m_{MH,}_{des}')
plot(x_m_MH_des, y_m_MH_des);
title('Normal distribution of m_{MH,}_{des}');
grid on;

subplot(4,2,3);
plot(x_T_abs, y_T_abs);
title('Normal distribution of T_{abs}');
grid on;

subplot(4,2,4);
plot(x_T_des, y_T_des);
title('Normal distribution of T_{des}');
grid on;

subplot(4,2,5);
plot(x_Peq_abs, y_Peq_abs);
title('Normal distribution of P_{eq,}_{abs}');
grid on;

subplot(4,2,6);
plot(x_Peq_des, y_Peq_des);
title('Normal distribution of P_{eq,}_{des}');
grid on;

subplot(4,2,7);
plot(x_r_abs, y_r_abs);
title('Normal distribution of r_{abs}');
grid on;

subplot(4,2,8);
plot(x_r_des, y_r_des);
title('Normal distribution of r_{des}');
grid on;

%% Other errors calculation (RMSE, MAE, ...)

% m_MH absorption
nan_m_MH_abs = isnan(m_MH_a_interp);
m_MH_a_interp(nan_m_MH_abs) = 0;
rmse_m_MH_a = sqrt(mean((m_MH_a_interp - m_MH_abs_paper).^2));
mae_m_MH_a = mean(abs(m_MH_a_interp - m_MH_abs_paper));
max_err_m_MH_a = max(abs(m_MH_a_interp - m_MH_abs_paper));
rel_err_m_MH_a = mean(abs(m_MH_a_interp - m_MH_abs_paper)) / mean(m_MH_abs_paper);
percent_err_m_MH_a = mean(abs(m_MH_a_interp - m_MH_abs_paper)) / mean(m_MH_abs_paper) * 100;
% m_MH desorption
nan_m_MH_des = isnan(m_MH_d_interp);
m_MH_d_interp(nan_m_MH_des) = 0;
rmse_m_MH_d= sqrt(mean((m_MH_d_interp - m_MH_des_paper).^2));
mae_m_MH_d = mean(abs(m_MH_d_interp - m_MH_des_paper));
max_err_m_MH_d = max(abs(m_MH_d_interp - m_MH_des_paper));
rel_err_m_MH_d = mean(abs(m_MH_d_interp - m_MH_des_paper)) / mean(m_MH_des_paper);
percent_err_m_MH_d = mean(abs(m_MH_d_interp - m_MH_des_paper)) / mean(m_MH_des_paper) * 100;


% T absorption
nan_T_abs = isnan(T_a_interp);
T_a_interp(nan_T_abs) = 0;
rmse_T_a = sqrt(mean((T_a_interp - T_abs_paper).^2));
mae_T_a = mean(abs(T_a_interp - T_abs_paper));
max_err_T_a = max(abs(T_a_interp - T_abs_paper));
rel_err_T_a = mean(abs(T_a_interp - T_abs_paper)) / mean(T_abs_paper);
percent_err_T_a = mean(abs(T_a_interp - T_abs_paper)) / mean(T_abs_paper) * 100;
% T desorption
nan_T_des = isnan(T_d_interp);
T_d_interp(nan_T_des) = 0;
rmse_T_d= sqrt(mean((T_d_interp - T_des_paper).^2));
mae_T_d = mean(abs(T_d_interp - T_des_paper));
max_err_T_d = max(abs(T_d_interp - T_des_paper));
rel_err_T_d = mean(abs(T_d_interp - T_des_paper)) / mean(T_des_paper);
percent_err_T_d = mean(abs(T_d_interp - T_des_paper)) / mean(T_des_paper) * 100;


% Peq absorption
nan_Peq_abs = isnan(Peq_a_interp);
Peq_a_interp(nan_Peq_abs) = 0;
rmse_Peq_a = sqrt(mean(((Peq_a_interp/1E+05) - Peq_abs_paper).^2));
mae_Peq_a = mean(abs((Peq_a_interp/1E+05) - Peq_abs_paper));
max_err_Peq_a = max(abs((Peq_a_interp/1E+05) - Peq_abs_paper));
rel_err_Peq_a = mean(abs((Peq_a_interp/1E+05) - Peq_abs_paper)) / mean(Peq_abs_paper);
percent_err_Peq_a = mean(abs((Peq_a_interp/1E+05) - Peq_abs_paper)) / mean(Peq_abs_paper) * 100;
% Peq desorption
nan_Peq_des = isnan(Peq_d_interp);
Peq_d_interp(nan_Peq_des) = 0;
rmse_Peq_d= sqrt(mean(((Peq_d_interp/1E+05) - Peq_des_paper).^2));
mae_Peq_d = mean(abs((Peq_d_interp/1E+05) - Peq_des_paper));
max_err_Peq_d = max(abs((Peq_d_interp/1E+05) - Peq_des_paper));
rel_err_Peq_d = mean(abs((Peq_d_interp/1E+05) - Peq_des_paper)) / mean(Peq_des_paper);
percent_err_Peq_d = mean(abs((Peq_d_interp/1E+05) - Peq_des_paper)) / mean(Peq_des_paper) * 100;


% r absorption
nan_r_abs = isnan(r_a_interp);
r_a_interp(nan_r_abs) = 0;
rmse_r_a = sqrt(mean((r_a_interp - r_abs_paper).^2));
mae_r_a = mean(abs(r_a_interp - r_abs_paper));
max_err_r_a = max(abs(r_a_interp - r_abs_paper));
rel_err_r_a = mean(abs(r_a_interp - r_abs_paper)) / mean(r_abs_paper);
percent_err_r_a = mean(abs(r_a_interp - r_abs_paper)) / mean(r_abs_paper) * 100;
% r desorption
nan_r_des = isnan(r_d_interp);
r_d_interp(nan_r_des) = 0;
rmse_r_d= sqrt(mean((r_d_interp - r_des_paper).^2));
mae_r_d = mean(abs(r_d_interp - r_des_paper));
max_err_r_d = max(abs(r_d_interp - r_des_paper));
rel_err_r_d = mean(abs(r_d_interp - r_des_paper)) / mean(r_des_paper);
percent_err_r_d = mean(abs(r_d_interp - r_des_paper)) / mean(r_des_paper) * 100;

% figure
% subplot(2,1,1);
% plot(ta, Peq_poly/1E+05);
% ylabel('Peq [bar]');
% xlabel('Time [s]');
% subplot(2,1,2);
% plot(H_M_a, Peq_poly);
% xlabel('H/M');
% ylabel('Peq [bar]');

%% Calcolo portata di acqua
Q = U*A*(T_a-Twa)/1000; %kW
Cp_water = 4.185; % kJ/kg K 
DeltaT_water = 5; % K

figure
plot(ta,Q);
ylabel('Q [kW]');
ylim([-5 35]);
xlabel('Time [s]');
title('Q_{out,reactor} => (Reactor ---> Water)');
xlim([-100 1800]);

Qmean = mean(Q);
m_w = 25/(Cp_water*DeltaT_water);

Q_des = abs(U*A*(T_d-Twd)/1000);

figure
plot(td,Q_des);
ylabel('Q [kW]');
ylim([-5 35]);
xlabel('Time [s]');
title('Q_{out,reactor} => (Reactor ---> Water)');
xlim([1700 3600]);

toc