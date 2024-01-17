clear all; clc;

%% ASSORBIMENTO
% Parametri assorbimento
Vg_abs = 0.0172;        % m3            comune
Cabs = 59.2;            % 1/s
Eabs = 21170;           % J/mol
R_abs = 8.314;          % J/ mol K      comune
P_abs = 5E+05;          % Pa
m_s_abs = 143;          % kg
M_H2_a = 2.016/1000;    % kg/mol        comune
SC_a = 3;               % molH2/molMH   comune
M_MH_a = 0.432;         % kg/mol        comune
Cp_H2_abs = 14300;      % J/kg K        comune
Cp_s_abs = 355;         % J/kg K        comune
Tin_abs = 290;          % K
A_abs = 5.4;            % m2            comune
U_abs = 243;            % W/m2 K        comune
Twa = 298;              % K
DHa = -30478;           % J/mol
DSa = -108;             % J/mol K
sl_abs = 0.13;          % -             comune
P0_abs = 1E+05;         % Pa            comune
m_H2_in = 2/3600;       % kg/s

%Condizioni iniziali
m_H2_0_abs = 0;         % kg
m_MH_0_abs = 0;         % kg
T_0_abs = Twa;          % K

y0_abs = [m_H2_0_abs; m_MH_0_abs; T_0_abs];

%Tempo di simulazione tspan
t_start_abs = 0;
t_end_abs = 1800;

tspan_abs = [t_start_abs, t_end_abs];

%Function handle
ode_fun_abs = @(tabs, yabs) Lumped_Model_caso1abs(tabs, yabs, m_H2_in);

%Risoluzione del sistema di equazioni
[tabs, yabs] = ode45(ode_fun_abs, tspan_abs, y0_abs);

%Estrazione delle variabili
m_H2_abs = yabs(:,1); % kg
m_MH_abs = yabs(:,2); % kg
T_abs = yabs(:,3);    % K

%Equazioni addizionali
Peq_abs = exp((DHa./(R_abs.*T_abs))-DSa./R_abs+sl_abs.*((m_MH_abs./m_s_abs)-0.5))*P0_abs;   % Pa
r_abs = Cabs.*exp(-Eabs./(R_abs.*T_abs)).*log(P_abs./Peq_abs).*(1-(m_MH_abs./m_s_abs));     % Reaction rate
Pg_abs = (m_H2_abs*R_abs.*T_abs./(Vg_abs.*(M_H2_a/1000)));                                  % Pa
wt_abs = (m_MH_abs./m_s_abs)*(M_H2_a*SC_a/M_MH_a)*100;                                      % [%]
H_M_abs = (wt_abs.*(M_MH_a/M_H2_a))./(100.-wt_abs);

%% DESORBIMENTO
% Parametri desorbimento
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

% Condizioni iniziali
m_H2_des_0 = 0;         % kg
m_MH_des_0 = m_s_des;   % kg
T_des_0 = Twd;          % K

y0_des = [m_H2_des_0; m_MH_des_0; T_des_0];

% Tempo di simulazione tspan
t_start_des = t_end_abs+1;
t_end_des = 3600;

tspan_des = [t_start_des, t_end_des];

% Function handle
ode_fun = @(t, y) Lumped_Model_caso1des(t, y);

% Risoluzione del sistema di equazioni
[t, y] = ode45(ode_fun, tspan_des, y0_des);

% Estrazione delle variabili
m_H2_des = y(:,1); % kg
m_MH_des = y(:,2); % kg
T_des = y(:,3);    % K

% Equazioni addizionali
Peq_des = exp(-(DHd./(R_des.*T_des))+DSd./R_des+sl_des.*((m_MH_des./m_s_des)-0.5))*P0_des;  % Pa
r_des = Cdes.*exp(-Edes./(R_des.*T_des)).*((P_des-Peq_des)./Peq_des).*(m_MH_des./m_s_des);  % Reaction rate
Pg_des = (m_H2_des*R_des.*T_des./(Vg_des.*(M_H2_des/1000)));                                % Pa
wt_des = (m_MH_des./m_s_des)*(M_H2_des*SC_des/M_MH_des)*100;                                % [%]
H_M_des = (wt_des.*(M_MH_des/M_H2_des))./(100.-wt_des);
 
%% Plot the results
figure(1)
subplot(4, 1, 1);
yyaxis left
plot(tabs, m_MH_abs);
ylabel('m_M_H [kg]');
ylim([-30 180]);
hold on
yyaxis right
plot(t, m_MH_des);
ylabel('m_M_H [kg]');
ylim([-30 180]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

subplot(4, 1, 2);
yyaxis left
plot(tabs, T_abs);
ylabel('T [K]');
ylim([290 360]);
hold on
yyaxis right
plot(t, T_des);
ylabel('T [K]');
ylim([290 360]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

subplot(4, 1, 3);
yyaxis left
plot(tabs, Peq_abs/1E+05);
ylabel('P_e_q [bar]');
ylim([0 15]);
hold on
yyaxis right
plot(t, Peq_des/1E+05);
ylabel('P_e_q [bar]');
ylim([0 15]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

subplot(4, 1, 4);
yyaxis left
plot(tabs, r_abs);
ylabel('r [1/s]');
ylim([-0.0030 0.0030]);
hold on
yyaxis right
plot(t, r_des);
ylabel('r [1/s]');
ylim([-0.0030 0.0030]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

figure(2)
subplot(3, 1, 1);
yyaxis left
plot(tabs, wt_abs);
ylabel('wt [%]');
ylim([0 2]);
hold on
yyaxis right
plot(t, wt_des);
ylabel('wt [%]');
ylim([0 2]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

subplot(3, 1, 2);
yyaxis left
plot(tabs, m_H2_abs);
ylabel('m_H_2 [kg]');
ylim([-5 5]);
hold on
yyaxis right
plot(t, m_H2_des);
ylabel('m_H_2 [kg]');
ylim([-5 5]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

subplot(3, 1, 3);
yyaxis left
plot(tabs, Pg_abs/1E+05);
ylabel('P_g [bar]');
%ylim([-10 0]);
hold on
yyaxis right
plot(t, Pg_des/1E+05);
ylabel('P_g [bar]');
%ylim([-5 5]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

figure(3)
yyaxis left
plot(tabs, H_M_abs);
ylabel('H/M');
ylim([0 5]);
hold on 
yyaxis right
plot(t, H_M_des);
ylabel('H/M');
ylim([0 5]);
hold off
xlabel('Time [s]');
xlim([0 3600]);
grid on

% 
% figure(3)
% plot(wt_des,Peq_des);
% xlabel('x')
% ylabel('P_e_q [Pa]')
% grid on 
% hold on
