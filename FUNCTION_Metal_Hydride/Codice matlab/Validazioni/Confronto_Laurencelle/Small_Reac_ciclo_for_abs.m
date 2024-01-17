clear all; clc;

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
m_H2_in = ((1:29) * 1E-09);    % kg/s
m_H2_in = [m_H2_in, (1:9) * 1E-08, (1:9) * 1E-07, 1E-06, 2E-06, 3E-06, 4E-06, 5E-06, 6E-06, 7E-06, 8E-06, 9E-06, 1E-05]';
Twa = 296;          % K

for i = 1:length(m_H2_in)
    % Condizioni iniziali
    m_H2_0 = 0;
    m_MH_0 = 0;
    T0 = Twa;
    
    y0 = [m_H2_0; m_MH_0; T0];
    
    % Tempo di simulazione 
    t_in = 0;
    t_fine = 900; 
    tspan = [t_in, t_fine];
    
    ode_fun = @(ta, ya) Abs_Laur_Small_Reactor(ta, ya, m_H2_in(i));
    
    [t, y] = ode45(ode_fun, tspan, y0);
    
    % Estrazione delle variabili
    m_H2_a = y(:, 1); % kg
    m_MH_a = y(:, 2); % kg
    T_a = y(:, 3);    % K
    
    % Equazioni addizionali
    Peq_abs = exp((DHa./(R.*T_a))-DSa./R+sl.*((m_MH_a./m_s)));     % Reaction rate
    Pg_abs = abs(m_H2_a.*R.*T_a./(Vg_a.*M_H2));                                      % Pa
    wt_abs = (m_MH_a./m_s)*(M_H2*SC/M_MH)*100;                                      % [%]
    H_M_abs = (wt_abs*(M_MH/M_H2))./(100-wt_abs);

    % Plot dei risultati
    figure(1)
    subplot(2, 1, 1);
    plot(t, wt_abs);
    xlabel('Time');
    ylabel('wt [%]');
    axis([0 900 0 1.5]);
    grid on
    hold on

    subplot(2, 1, 2);
    plot(t, T_a);
    xlabel('Time');
    ylabel('T [K]');
    axis([0 900 270 340]);
    grid on
    hold on

    % Aggiunta della legenda
    legendInfo{i} = ['m_{H2,}_{in} = ', num2str(m_H2_in(i))];
end

% Aggiunta della legenda finale
figure(1)
subplot(2, 1, 1);
legend(legendInfo);
subplot(2, 1, 2);
legend(legendInfo);
