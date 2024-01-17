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
Twa      = 280:2:294;  % Inlet Cooling water temperature [K]
Twd      = 353;     % Inlet Heating water temperature    [K]

%% ABSORPTION

%Absorption parameters
Ca    = 59.2;   % Kinetic constant                   [1/s]
Ea    = 21170;  % Absorption activation energy       [J/mol]
Pa    = 5E+05;  % Absorption pressure                [Pa]
Tin_a = 290;    % H2 input temperature               [K]
DHa   = -30478; % Absorption enthalpy                [J/mol]
DSa   = -108;   % Absorption entropy                 [J/mol K]
m_w   = 0.5:1:10.5;  % Portata d'acqua di progetto   [kg/s]
Cp_w  = 4.185;       % Specific heat water           [kJ/kg K] 

t = (0:1:1800)'; % Tempo di simulazione 

% Figure 1: Results at varying m_w
figure(1);
subplot(2, 1, 1);
hold on;
grid on;
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylabel('T_a [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
title('System Temperature at Varying m_w', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
legendInfo = {};

subplot(2, 1, 2);
hold on;
grid on;
xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
ylabel('T_{w,out} [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
title('Outlet Water Temperature at Varying m_w', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
legendInfo = {};

for j = 1:numel(m_w)
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
    T0_a    = Twa(1)+15;  % System's temperature            [K]
    Tw_out0 = Twa(1)+5;

    % Initial conditions for the first timestep
    m_H2_a(1) = m_H20_a;
    m_MH_a(1) = m_MH0_a;
    T_a(1)    = T0_a;
    Tw_out(1) = Tw_out0;
    DTb(1) = T_a(1) - Twa(1);
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
        ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in, Twa(1), DTml(i));

        % Differential equation solver
        [ta, ya] = ode23s(ode_fun_a, tspan, y0);

        % Store solutions for current timestep 
        m_H2_a(i+1) = ya(end, 1);
        m_MH_a(i+1) = ya(end, 2);
        T_a(i+1)    = ya(end, 3);

        % Calculate outlet water temperature
        DTb(i+1) = T_a(i) - Twa(1);
        DTa(i+1) = T_a(i) - Tw_out(i);
        DTml(i+1) = (DTb(i+1) - DTa(i+1)) / log(DTb(i+1) / DTa(i+1));
        Tw_out(i+1) = Twa(1) + (DTml(i+1) * U * A) / (m_w(j) * Cp_w * 1000);
    end

    subplot(2, 1, 1);
    plot(t, T_a);
    legendInfo{j} = sprintf('m_w = %.1f kg/s', m_w(j));

    subplot(2, 1, 2);
    plot(t, Tw_out);
    legendInfo{j} = sprintf('m_w = %.1f kg/s', m_w(j));
end

subplot(2, 1, 1);
legend(legendInfo, 'Location', 'northeast');

subplot(2, 1, 2);
legend(legendInfo, 'Location', 'northeast');

% % Figures at varying Twa
% for k = 1:numel(Twa)
%     figure(k+1);
%     subplot(2, 1, 1);
%     hold on;
%     grid on;
%     xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     ylabel('T_a [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     title(sprintf('System Temperature (Twa = %.1f K)', Twa(k)), 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     legendInfo = {};
% 
%     subplot(2, 1, 2);
%     hold on;
%     grid on;
%     xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     ylabel('T_{w,out} [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     title(sprintf('Outlet Water Temperature (Twa = %.1f K)', Twa(k)), 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     legendInfo = {};
% 
%     for j = 1:numel(m_w)
%         % Inizializzazione variabili 
%         m_H2_a = zeros(size(t));
%         m_MH_a = zeros(size(t));
%         T_a    = zeros(size(t));
%         Tw_out = zeros(size(t));
%         DTb    = zeros(size(t));
%         DTa    = zeros(size(t));
%         DTml   = zeros(size(t));
% 
%         %Initial conditions
%         m_H20_a = 0;       % Absorbed hydrogen mass          [kg]
%         m_MH0_a = 0;       % Metal hydride mass              [kg]
%         T0_a    = Twa(k)+15;  % System's temperature            [K]
%         Tw_out0 = Twa(k)+5;
% 
%         % Initial conditions for the first timestep
%         m_H2_a(1) = m_H20_a;
%         m_MH_a(1) = m_MH0_a;
%         T_a(1)    = T0_a;
%         Tw_out(1) = Tw_out0;
%         DTb(1) = T_a(1) - Twa(k);
%         DTa(1) = T_a(1) - Tw_out(1);
%         DTml(1) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));
% 
%         for i = 1:numel(t)-1
%             if i == 1
%                 y0 = [m_H2_a(1); m_MH_a(1); T_a(1)];
%                 DTb(i) = DTb(1);
%                 DTa(i) = DTa(1);
%                 DTml(i) = (DTb(1)-DTa(1))/(log(DTb(1)/DTa(1)));
%             else
%                 y0 = [m_H2_a(i); m_MH_a(i); T_a(i)];
%             end
% 
%             % Define tspan for ODE resolution for current timestep
%             tspan = [t(i) t(i+1)];
% 
%             % Function handle
%             ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in, Twa(k), DTml(i));
% 
%             % Differential equation solver
%             [ta, ya] = ode23s(ode_fun_a, tspan, y0);
% 
%             % Store solutions for current timestep 
%             m_H2_a(i+1) = ya(end, 1);
%             m_MH_a(i+1) = ya(end, 2);
%             T_a(i+1)    = ya(end, 3);
% 
%             % Calculate outlet water temperature
%             DTb(i+1) = T_a(i) - Twa(k);
%             DTa(i+1) = T_a(i) - Tw_out(i);
%             DTml(i+1) = (DTb(i+1) - DTa(i+1)) / log(DTb(i+1) / DTa(i+1));
%             Tw_out(i+1) = Twa(k) + (DTml(i+1) * U * A) / (m_w(j) * Cp_w * 1000);
%         end
% 
%         subplot(2, 1, 1);
%         plot(t, T_a);
%         legendInfo{j} = sprintf('m_w = %.1f kg/s', m_w(j));
% 
%         subplot(2, 1, 2);
%         plot(t, Tw_out);
%         legendInfo{j} = sprintf('m_w = %.1f kg/s', m_w(j));
%     end
% 
%     subplot(2, 1, 1);
%     legend(legendInfo, 'Location', 'northeast');
% 
%     subplot(2, 1, 2);
%     legend(legendInfo, 'Location', 'northeast');
% end

