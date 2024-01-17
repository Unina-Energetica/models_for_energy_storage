clear all; clc;

%% MODEL DESCRIPTION
% In this model I validate the equations by comparing the results I
% obtained, with the results of the SMALL REACTOR defined in the paper 

%% ASSORBIMENTO
% Parametri
Vg_a = 4.4E-07;     % m3
Ca = 59.2;          % 1/s
Ea = 21170;         % J/mol
R = 8.314;          % J/ mol K
Pa = 6E+05;         % Pa
m_s = 1/1000;       % kg
M_H2 = 2.016/1000;  % kg/mol
SC = 2.76;          % molH2/molMH
M_MH = 0.432;       % kg/mol
Cp_H2 = 14300;      % J/kg K
Cp_s = 355;         % J/kg K 
Tin_a = 290;        % K
A = 5E-04;          % m2
U = 80;             % W/m2 K
DHa = -30478;       % J/mol
DSa = -108;         % J/mol K
sl = 0.13;          % -
P0 = 1E+05;         % Pa
Twa = 296;          % K
m_H2_in = 1E-09;    % kg/s

% Condizioni iniziali
    m_H2_0 = 0;
    m_MH_0 = 0;
    T0 = Twa;
    
    y0 = [m_H2_0; m_MH_0; T0];
    
    % Tempo di simulazione 
    t_in = 0;
    t_fine = 900; 
    tspan = [t_in, t_fine];
    
    ode_fun = @(ta, ya) Abs_Laur_Small_Reactor(ta, ya, m_H2_in);
    
    [t, y] = ode23s(ode_fun, tspan, y0);
    
    % Estrazione delle variabili
    m_H2_a = y(:, 1); % kg
    m_MH_a = y(:, 2); % kg
    T_a = y(:, 3);    % K
    
    % Equazioni addizionali
    Peq_abs = exp((DHa./(R.*T_a))-DSa./R+sl.*((m_MH_a./m_s)-0.5))*P0;     % Reaction rate
    Pg_abs = abs(m_H2_a.*R.*T_a./(Vg_a.*M_H2));                                      % Pa
    wt_abs = (m_MH_a./m_s)*(M_H2*SC/M_MH)*100;                                      % [%]
    H_M_abs = (wt_abs*(M_MH/M_H2))./(100-wt_abs);



%% DESORBIMENTO 
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
m_H2_out = 1E-09;    % kg/s

%Condizioni iniziali 
m_H20_d = 0;
m_MH0_d = m_s;
T0_d = Twd;

y0_d = [m_H20_d; m_MH0_d; T0_d];


  %Tempo di simulazione
    t_start_d = t_fine;
    t_end_d = 1800;

    t_span_d = [t_start_d, t_end_d];

    %Function handle 
    ode_fun_d = @(td, yd) Des_Laur_Small_Reactor(td, yd, m_H2_out);

    %Risoluzione del sistema di equazioni
    [td, yd] = ode23s(ode_fun_d,t_span_d,y0_d);

    %Estrazione delle variabili
    m_H2_d = yd(:,1);   % [kg]
    m_MH_d = yd(:,2);   % [kg]
    T_d = yd(:,3);      % [K]

    % Equazioni addizionali
    Peq_d = exp(-(DHd./(R.*T_d))+DSd./R+sl.*((m_MH_d./m_s)-0.5))*P0; %Pressione equilibrio desorbimento [Pa]
    r_d = Cd.*exp(-Ed./(R.*T_d)).*((Pd-Peq_d)./Peq_d).*(m_MH_d./m_s);  %Reaction rate                     [1/s]
    Pg_d = (m_H2_d*R.*T_d./(Vg_a.*(M_H2/1000)));                         %Percentuale assorbita             [Pa]
    wt_d = (m_MH_d./m_s)*(M_H2*SC/M_MH)*100;                           %Rapporto concentrazioni           [-]
    H_M_d = (wt_d.*(M_MH/M_H2))./(100.-wt_d);

%% RISULTATI 
% % Plot dei risultati
% figure(1)
% subplot(2, 1, 1);
% yyaxis left
% plot(t, wt_abs);
% ylabel('wt [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
% ylim([0 2]);
% hold on
% yyaxis right
% plot(td, wt_d);
% ylabel('wt [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
% ylim([0 2]);
% hold off
% xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
% xticks(0:100:1800);
% grid on
% 
% subplot(2,1,2);
% yyaxis left
% plot(t, T_a-273);
% ylabel('T [°C]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
% ylim([-20 60]);
% hold on
% yyaxis right
% plot(td, T_d-273);
% ylabel('T [°C]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
% ylim([-20 60]);
% hold off
% xlabel('Time [s]');
% xticks(0:100:1800);
% grid on

time = cat(1, t, td);
Pressure = cat(1, Peq_abs, Peq_d);

figure(2)
plot(time, Pressure/1E+05);
ylabel('P [bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylim([0 7]);
xlabel('Time [s]');
xticks(0:100:1800);
grid on


%% IMPORTAZIONE ED INTERPOLAZIONE RISULTATI LAURENCELLE

%Temperatura 
ta_Laur = xlsread('T_Laur.xlsx','Foglio1','A1:A130');
Ta_Laur = xlsread('T_Laur.xlsx', 'Foglio1', 'B1:B130');
td_Laur = xlsread('T_Laur_Des.xlsx','Foglio1','A1:A183');
Td_Laur = xlsread('T_Laur_Des.xlsx','Foglio1','B1:B183');

%Interpolazione temperatura
Ta_interp = interp1(t, T_a, ta_Laur);
Td_interp = interp1(td, T_d, td_Laur);

% wt modello
ta_wt_Laur = xlsread('wt_abs_Laur.xlsx','Foglio1','A1:A198');
wt_a_Laur = xlsread('wt_abs_Laur.xlsx','Foglio1','B1:B198');
td_wt_Laur = xlsread('wt_des_Laure.xlsx','Foglio1','A1:A169');
wt_d_Laur = xlsread('wt_des_Laure.xlsx','Foglio1','B1:B169');

% wt sperimentale
wt_1_exp = readtable("wt_a_Laur_exper.xlsx");
ta_wt_exp = wt_1_exp.Var1;
wt_a_exp = wt_1_exp.Var2;
wt_2_exp = readtable("wt_d_Laur_exper.xlsx");
td_wt_exp = wt_2_exp.Var1;
wt_d_exp = wt_2_exp.Var2;

% Interpolazione wt
wt_a_interp = interp1(t, wt_abs, ta_wt_Laur);
wt_d_interp = interp1(td, wt_d, td_wt_Laur);
wt_a_exp_interp = interp1(ta_wt_exp, wt_a_exp, ta_wt_Laur);
wt_d_exp_interp = interp1(td_wt_exp, wt_d_exp, td_wt_Laur);


%% UNIONE DI TUTTI I VETTORI CREATI E PLOT

% Unione vettori
tempo_T = cat(1, ta_Laur, td_Laur);
tempo_wt = cat(1, ta_wt_Laur, td_wt_Laur);
Temp_Modello = cat(1, Ta_interp, Td_interp);
Temp_Laurencelle = cat(1, Ta_Laur, Td_Laur);
wt_Modello = cat(1, wt_a_interp, wt_d_interp);
wt_Laurencelle = cat(1, wt_a_Laur, wt_d_Laur);
wt_Experimental = cat(1, wt_a_exp_interp, wt_d_exp_interp);

figure(3)
subplot(2,1,1);
yyaxis left
plot(tempo_T, Temp_Modello-273, 'b');
ylabel('T_{Antonio} [°C]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([-20 60]);
yyaxis right 
plot(tempo_T, Temp_Laurencelle-273, 'r');
ylabel('T_{Laurencelle} [°C]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([-20 60]);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
xticks(0:100:1800);
grid on
legend('T_{Antonio}', 'T_{Laurencelle}', 'location', 'northeast', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);

subplot(2,1,2);
yyaxis left
plot(tempo_wt, wt_Modello, 'b');
ylabel('wt [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([0 1.4]);
yyaxis right 
plot(tempo_wt, wt_Laurencelle, 'r', tempo_wt, wt_Experimental, 'c-');
ylabel('wt [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([0 1.4]);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
xticks(0:100:1800);
grid on
legend('wt_{Antonio}', 'wt_{Laurencelle}', 'wt_{Experimental}', 'location','northeast', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);


% figure(2)
% subplot(2,1,1)
% yyaxis left
% plot(ta_Laur, Ta_Laur-273, 'b-', td_Laur, Td_Laur-273, 'b--');
% ylabel('T_{Laurencelle} [K]');
% ylim([-20 60]);
% yyaxis right
% plot(ta_Laur, Ta_interp-273, 'r-', td_Laur, Td_interp-273, 'r--');
% ylabel('T_{Antonio} [K]');
% ylim([-20 60]);
% xlabel('Time [s]');
% xlim([0 1800]);
% grid on
% legend('T_a_{,Laurencelle}', 'T_d_{,Laurencelle}', 'T_a_{,Antonio}', 'T_d_{,Antonio}','location','northeast');
% 
% subplot(2,1,2);
% yyaxis left
% plot(ta_wt_Laur, wt_a_Laur, 'b-', td_wt_Laur, wt_d_Laur, 'b--', ta_wt_Laur, wt_a_exp_interp, 'c-', td_wt_Laur, wt_d_exp_interp, 'c--');
% ylabel('wt [%]');
% ylim([0 1.4]);
% yyaxis right
% plot(ta_wt_Laur, wt_a_interp, 'r-', td_wt_Laur, wt_d_interp, 'r--');
% ylabel('wt [%]');
% ylim([0 1.4]);
% xlabel('Time [s]');
% xlim([0 1800]);
% grid on
% legend('wt_{abs,Laurencelle}', 'wt_{des,Laurencelle}', 'wt_{abs,Experimental}', 'wt_{des,Experimental}', 'wt_{abs,Antonio}','wt_{des,Antonio}','location','northeast');


%% CALCOLO ERRORE PERCENTUALE E PLOT

% T
err_T_Laur = abs(((Temp_Modello - Temp_Laurencelle) ./ Temp_Modello) * 100);

% wt
err_wt_Laur = abs(((wt_Modello - wt_Laurencelle) ./ wt_Modello) * 100);
err_wt_Exp = abs(((wt_Modello - wt_Experimental) ./ wt_Modello) * 100);

figure(4)
subplot(2, 1, 1);
plot(tempo_T, err_T_Laur);
ylabel('errore T', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([-5 100]);
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
xticks(0:100:1800);
grid on

subplot(2,1,2);
yyaxis left 
plot(tempo_wt,err_wt_Laur);
ylabel('errore_{wt,Laurencelle}', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([-5 100]);
hold on
yyaxis right
plot(tempo_wt, err_wt_Exp);
ylabel('errore_{wt,Experimental}', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);
ylim([-5 100]);
hold off
xlabel('Time [s]');
xticks(0:100:1800);
grid on
legend('Errore_{Laurencelle}','Errore_{Experimental}', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',16);