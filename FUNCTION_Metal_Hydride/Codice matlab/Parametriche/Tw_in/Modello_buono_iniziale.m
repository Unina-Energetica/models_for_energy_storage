clear all; clc;

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
Twa      = 298-15;  % Inlet Cooling water temperature    [K]
Twd      = 353;     % Inlet Heating water temperature    [K]

%% Parametrizzazione m_w

%Absorption parameters
Ca    = 59.2;   % Kinetic constant                   [1/s]
Ea    = 21170;  % Absorption activation energy       [J/mol]
Pa    = 5E+05;  % Absorption pressure                [Pa]
Tin_a = 290;    % H2 input temperature               [K]
DHa   = -30478; % Absorption enthalpy                [J/mol]
DSa   = -108;   % Absorption entropy                 [J/mol K]
m_w   = 1.1947;       % Portata d'acqua di progetto   [kg/s]  0.5:1:10.5
Cp_w  = 4.185;       % Specific heat water           [kJ/kg K] 

t = (0:1:1800)'; % Tempo di simulazione 

%for j= 1:numel(m_w)

% Inizializzazione variabili 
m_H2_a = zeros(size(t));
m_MH_a = zeros(size(t));
T_a    = zeros(size(t));
Tw_out = zeros(size(t));
DTb    = zeros(size(t));
DTa    = zeros(size(t));
DTml   = zeros(size(t));

%Initial conditions
m_H20_a = 0;       % Absorbed hydrogen mass          [kg]
m_MH0_a = 0;       % Metal hydride mass              [kg]
T0_a    = Twa+15;  % System's temperature            [K]
Tw_out0 = Twa+5;

% Initial conditions for the first timestep
m_H2_a(1) = m_H20_a;
m_MH_a(1) = m_MH0_a;
T_a(1)    = T0_a;
Tw_out(1) = Tw_out0;
DTb(1) = T_a(1) - Twa;
DTa(1) = T_a(1) - Tw_out(1);
DTml(1) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));


for i = 1:numel(t)-1

    if i == 1
        y0 = [m_H2_a(1); m_MH_a(1); T_a(1)];
        DTb(i) = DTb(1);
        DTa(i) = DTa(1);
        DTml(i) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));
    else
        y0 = [m_H2_a(i); m_MH_a(i); T_a(i)];
    end

    % Define tspan for ODE resolution for current timestep
    tspan = [t(i) t(i+1)];

    % Function handle
    ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in, Twa, DTml(i));
    
    % Differential equation solver
    [ta, ya] = ode23s(ode_fun_a, tspan, y0);
    
    % Store solutions for current timestep 
    m_H2_a(i+1) = ya(end, 1);
    m_MH_a(i+1) = ya(end, 2);
    T_a(i+1)    = ya(end, 3);

    % Calculate outlet water temperature
    Tw_out(i+1) = Twa + (DTml(i) * U * A) / (m_w * Cp_w * 1000);
    DTb(i+1) = T_a(i+1) - Twa;
    DTa(i+1) = T_a(i+1) - Tw_out(i+1);
    DTml(i+1) = (DTb(i+1) - DTa(i+1)) / log(DTb(i+1) / DTa(i+1));
    % Tw_out(i+1) = Twa + (DTml(i+1) * U * A) / (m_w(j) * Cp_w * 1000);
end

figure(1)
    subplot(2, 1, 1);
    plot(t, T_a);
    ylabel('T_a [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legendInfo = cellstr(num2str(m_w', 'm_{w,}_{in} = %.1f'));
    legend(legendInfo, 'Location', 'northeast');
    title('System Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

    subplot(2, 1, 2);
    plot(t, Tw_out);
    ylabel('T_{w,out} [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legendInfo = cellstr(num2str(m_w', 'm_{w,}_{in} = %.1f'));
    legend(legendInfo, 'Location', 'northeast');
    title('Outlet Water Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

%end

%% Parametrizzazione Tw_in
Twa      = 280:2:298;  % Inlet Cooling water temperature    [K]
m_w_pp = 1.1947;

for k= 1:numel(Twa)

% Inizializzazione variabili 
m_H2_a = zeros(size(t));
m_MH_a = zeros(size(t));
T_a    = zeros(size(t));
Tw_out = zeros(size(t));
DTb    = zeros(size(t));
DTa    = zeros(size(t));
DTml   = zeros(size(t));

%Initial conditions
m_H20_a = 0;       % Absorbed hydrogen mass          [kg]
m_MH0_a = 0;       % Metal hydride mass              [kg]
T0_a    = Twa(1)+3;  % System's temperature            [K]
Tw_out0 = Twa(1)+5;

% Initial conditions for the first timestep
m_H2_a(1) = m_H20_a;
m_MH_a(1) = m_MH0_a;
T_a(1)    = T0_a;
Tw_out(1) = Tw_out0;
DTb(1) = T_a(1) - (Twa(1)+3);
DTa(1) = T_a(1) - Tw_out(1);
DTml(1) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));


for i = 1:numel(t)-1
    if i == 1
        y0 = [m_H2_a(1); m_MH_a(1); T_a(1)];
        DTb(i) = DTb(1);
        DTa(i) = DTa(1);
        DTml(i) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));
    else
        y0 = [m_H2_a(i); m_MH_a(i); T_a(i)];
    end

    % Define tspan for ODE resolution for current timestep
    tspan = [t(i) t(i+1)];

    % Function handle
    ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in, Twa, DTml(i));

    % Differential equation solver
    [ta, ya] = ode23s(ode_fun_a, tspan, y0);

    % Store solutions for current timestep 
    m_H2_a(i+1) = ya(end, 1);
    m_MH_a(i+1) = ya(end, 2);
    T_a(i+1)    = ya(end, 3);

    % Calculate outlet water temperature
    DTb(i+1) = T_a(i) - Twa(k);
    DTa(i+1) = T_a(i) - Tw_out(i);
    DTml(i+1) = (DTb(i+1) - DTa(i+1)) / log(DTb(i+1) / DTa(i+1));
    Tw_out(i+1) = Twa(k) + (DTml(i+1) * U * A) / (m_w_pp * Cp_w * 1000);

end

figure(2)
    subplot(2, 1, 1);
    plot(t, T_a);
    ylabel('T_a [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legendInfo = cellstr(num2str(Twa', 'T_{w,}_{in} = %.1f'));
    legend(legendInfo, 'Location', 'northeast');
    title('System Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

    subplot(2, 1, 2);
    plot(t, Tw_out);
    ylabel('T_{w,out} [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legendInfo = cellstr(num2str(Twa', 'T_{w,}_{in} = %.1f'));
    legend(legendInfo, 'Location', 'northeast');
    title('Outlet Water Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

    figure(3)
    plot(t, m_MH_a);
    ylabel('m_{MH}_a [kg]','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legendInfo = cellstr(num2str(Twa', 'T_{w,}_{in} = %.1f'));
    legend(legendInfo, 'Location', 'northeast');
    title('Metal Hydride Mass', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on


end

% %% Plot per confronto con risultati vecchi
% Twa = 298;
% m_w = 1.1947;
% 
% for i = 1:numel(t)-1
%     if i == 1
%         y0 = [m_H2_a(1); m_MH_a(1); T_a(1)];
%         DTb(i) = DTb(1);
%         DTa(i) = DTa(1);
%         DTml(i) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));
%     else
%         y0 = [m_H2_a(i); m_MH_a(i); T_a(i)];
%     end
% 
%     % Define tspan for ODE resolution for current timestep
%     tspan = [t(i) t(i+1)];
% 
%     % Function handle
%     ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in, Twa, DTml(i));
% 
%     % Differential equation solver
%     [ta, ya] = ode23s(ode_fun_a, tspan, y0);
% 
%     % Store solutions for current timestep 
%     m_H2_a(i+1) = ya(end, 1);
%     m_MH_a(i+1) = ya(end, 2);
%     T_a(i+1)    = ya(end, 3);
% 
%     % Calculate outlet water temperature
%     DTb(i+1) = T_a(i) - Twa;
%     DTa(i+1) = T_a(i) - Tw_out(i);
%     DTml(i+1) = (DTb(i+1) - DTa(i+1)) / log(DTb(i+1) / DTa(i+1));
%     Tw_out(i+1) = Twa + (DTml(i+1) * U * A) / (m_w * Cp_w * 1000);
% 
% end
% 
% figure(300)
% plot(t, T_a);
% ylabel('T_a [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
% ylim([290 360]);
% title('System Temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
% xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
% hold on
% grid on


