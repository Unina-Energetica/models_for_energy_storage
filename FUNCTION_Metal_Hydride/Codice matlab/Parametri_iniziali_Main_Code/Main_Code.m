clear all; clc;
tic
%% Equations to be used
% Calculation of m_s based on the quantity of m_H2 to be absorbed
% m_s = (m_H2/wt_max) * 100  [kg]    con wt_max = 1.39
% 
% Calculation of hydrogen flowrate m_H2_in based on the quantity absorbed
% m_H2_in = (2 * m_H2) / 3600  [kg/s]

%% Parameters absorption and desorption
rho_MH = 8310;        % Density of LaNi5                   [kg/m3]
Vg     = 0.0172;      % Gas Volume                         [m3]
R      = 8.314;       % Gas constant                       [J/mol K]
M_H2   = 2.016/1000;  % Molar mass H2                      [kg/mol]
M_MH   = 0.432;       % Molar mass   MH                    [kg/mol]
Cp_H2  = 14300;       % Specific heat H2                   [J/kg K]
Cp_s   = 355;         % Specific heat solid                [J/kg K]
A      = 5.4;         % Area "heat interchange"            [m2]
U      = 243;         % Overall heat transfer coefficient  [W/m2 K]
sl     = 0.13;        % Slope coefficient                  [-]
P0     = 1E+05;       % Reference pressure                 [Pa]
m_s    = 143;         % Solid mass                         [kg]
SC     = 3;           % Stoichiometric coefficient         [-]
Cp_w   = 4.185;       % Specific heat water                [kJ/kg K] 

%% Simulation time for absorption and desorption
t_abs  = (0:1:1800)';                        % [s]
t_des  = (t_abs(end)+1:1:2*t_abs(end)+1)';   % [s]

%% Input parameters
m_H2_in   = ones(numel(t_abs),1)*(4/3600);   % Input hydrogen flowrate               [kg/s]
m_H2_out  = ones(numel(t_des),1)*(4/3600);   % Output hydrogen flowrate              [kg/s]
Tw_in_abs = ones(numel(t_abs),1) * 286;      % Inlet cooling water temperature       [k]
Tw_in_des = ones(numel(t_des),1) * 353;      % Inlet heating water temperature fixed [K]
m_w_abs   = ones(numel(t_abs),1) * 1.1947;   % Project water flowrate for absor      [kg/s]
m_w_des   = ones(numel(t_des),1) * 1.1947;   % Project water for desorption          [kg/s]

%% Absorption
%Absorption parameters
C_abs      = 59.2;   % Kinetic constant                   [1/s]
E_abs      = 21170;  % Absorption activation energy       [J/mol]
P_abs      = 8E+05;  % Absorption pressure                [Pa]
T_H2_in    = 290;    % H2 input temperature               [K]
DH_abs     = -30478; % Absorption enthalpy                [J/mol]
DS_abs     = -108;   % Absorption entropy                 [J/mol K]
eps_abs    = 10E-10; % Control value for energy balance   [-]

% Variable inizialization 
m_H2_abs   = zeros(size(t_abs));  % Absorbed hydrogen mass           [kg]
m_MH_abs   = zeros(size(t_abs));  % Metal hydride mass formed        [kg]
T_abs      = zeros(size(t_abs));  % System's temperature             [K]
Tw_out_abs = zeros(size(t_abs));  % Outlet water temperature         [K]
DTb_abs    = zeros(size(t_abs));  % 
DTa_abs    = zeros(size(t_abs));  % 
DTml_abs   = zeros(size(t_abs));  % 
Peq_abs    = zeros(size(t_abs));  % Equilibrium pressure             [Pa]
r_abs      = zeros(size(t_abs));  % Reaction rate                    [1/s]

% Initial conditions for the first timestep
m_H2_abs(1)   = 0;                        % Absorbed hydrogen mass           [kg]
m_MH_abs(1)   = 0;                        % Metal hydride mass               [kg]
T_abs(1)      = 298;                      % System's temperature             [K]
Tw_out_abs(1) = Tw_in_abs(1)+10;          % Initial outlet water temperature [K] 
DTb_abs(1)    = T_abs(1) - Tw_in_abs(1); 
DTa_abs(1)    = T_abs(1) - Tw_out_abs(1);  
DTml_abs(1)   = (DTb_abs(1)-DTa_abs(1))/(log(DTb_abs(1)/DTa_abs(1)));

for i = 1:numel(t_abs)-1    

    res_abs = 1000;         % Initial value for the residual

    while res_abs > eps_abs % We repeat the iteration if the energy balance is not verified

    if i == 1               %  "if" used to assign the initial condition based on the "i" index
        y0_a        = [m_H2_abs(1); m_MH_abs(1); T_abs(1)]; % Initial conditions for timestep = 1
        DTb_abs(i)  = DTb_abs(1);
        DTa_abs(i)  = DTa_abs(1);
        DTml_abs(i) = (DTb_abs(1)-DTa_abs(1))/(log(DTb_abs(1)/DTa_abs(1)));
    else
        y0_a        = [m_H2_abs(i); m_MH_abs(i); T_abs(i)]; % Initial conditions for timestep > 1
    end

    % Define tspan for ODE resolution for current timestep
    tspan = [t_abs(i) t_abs(i+1)];

    % Function handle
    ode_fun_a = @(ta, ya) Absorption(ta, ya, m_H2_in(i), Tw_in_abs(i), DTml_abs(i), m_s, A, U, P_abs, C_abs, E_abs, R, M_H2, SC, M_MH, Cp_H2, Cp_s, T_H2_in, DH_abs, DS_abs, sl, P0);

    % Differential equation solver
    [ta, ya] = ode23s(ode_fun_a, tspan, y0_a);

    % Store solutions for current timestep 
    m_H2_abs(i+1) = ya(end, 1);
    m_MH_abs(i+1) = ya(end, 2);
    T_abs(i+1)    = ya(end, 3);

    % Calculation of other parameters
    Peq_abs(i+1) = exp((DH_abs/(R*T_abs(i)))-DS_abs/R+sl*((m_MH_abs(i)/m_s)-0.5))*P0;           % Absorption equilibrium pressure [Pa]
    r_abs(i+1)   = C_abs*exp(-E_abs/(R*T_abs(i)))*log(P_abs/Peq_abs(i))*(1-(m_MH_abs(i)/m_s));  % Reaction rate                   [1/s]

    % Calculate outlet water temperature
    DTb_abs(i+1)    = T_abs(i) - Tw_in_abs(i);                                                  % Calculation of the new DTb               [K]
    DTa_abs(i+1)    = T_abs(i) - Tw_out_abs(i);                                                 % Calculation of the new DTa               [K]
    DTml_abs(i+1)   = (DTb_abs(i+1) - DTa_abs(i+1)) / log(DTb_abs(i+1) / DTa_abs(i+1));         % Calculation of the new DTml              [K]
    Tw_out_abs(i+1)   = Tw_in_abs(i) + (DTml_abs(i+1) * U * A) / (m_w_abs(i) * Cp_w * 1000);    % Calculation of the new Tw_out            [K]
    Q_abs(i+1)      = (m_w_abs(i)*Cp_w*(Tw_out_abs(i+1)-Tw_in_abs(i)));                         % Calculation of the heat (MH --> water)   [kW]

    res_abs         = DTml_abs(i+1)*U*A/1000 - (Tw_out_abs(i+1)-Tw_in_abs(i))*m_w_abs(i)*Cp_w;  % Energy balance equation to be verified 
    end
end

Qabs = Q_abs'; % Transposition of the row vector in a column vector 

%% Desorption
% Desorption parameters

C_des      = 9.6;      % Kinetic constant              [1/s]
E_des      = 16420;    % Desorption activation energy  [J/mol]
P_des      = 5E+05;    % Desorption pressure           [Pa]
DH_des     = 30800;    % Desorption enthalpy           [J/mol]
DS_des     = 108;      % Desorption entropy            [J/mol K]
eps_des    = 10E-10;   % Control value for energy balance

% Variable inizialization
m_H2_des   = zeros(size(t_des)); % Desorbed hydrogen mass   [kg]
m_MH_des   = zeros(size(t_des)); % Metal hydride mass left  [kg]
T_des      = zeros(size(t_des)); % System's temperature     [K]
Tw_out_des = zeros(size(t_des)); % Outlet water temperature
DTb_des    = zeros(size(t_des)); %
DTa_des    = zeros(size(t_des)); %
DTml_des   = zeros(size(t_des)); %
Peq_des    = zeros(size(t_des)); % Equilibrium pressure     [Pa]
r_des      = zeros(size(t_des)); % Reaction rate            [1/s]

% Initial conditions for the first timestep
m_H2_des(1)   = 0;                     % Desorbed hydrogen mass            [kg]
m_MH_des(1)   = m_MH_abs(end);         % Metal hydride mass left           [kg]
T_des(1)      = T_abs(end);            % System's temperature              [K]
Tw_out_des(1) = Tw_in_des(1)-10;       % Initial outlet water temperature  [K]
DTb_des(1)  = Tw_in_des(1) - T_des(1);
DTa_des(1)  = Tw_out_des(1) - T_des(1);
DTml_des(1) = (DTb_des(1)-DTa_des(1))/(log(DTb_des(1)/DTa_des(1)));
Peq_des(1)    = Peq_abs(end);          % Initial equilibrium pressure      [Pa]


for i = 1:numel(t_des)-1

    res_des = 1000;         % Initial value for the residual

    while res_des > eps_des % We repeat the iteration if the energy balance is not verified

        if i == 1           %  "if" used to assign the initial condition based on the "i" index
            y0_d        = [m_H2_des(1); m_MH_des(1); T_des(1)]; % Initial conditions for timestep = 1
            DTb_des(i)  = DTb_des(1);
            DTa_des(i)  = DTa_des(1);
            DTml_des(i) = DTml_des(1);
        else
            y0_d = [m_H2_des(i); m_MH_des(i); T_des(i)];        % Initial conditions for timestep > 1
        end

    % Define tspan for ODE resolution for current timestep
    tspan = [t_des(i) t_des(i+1)];

    % Function handle
    ode_fun_d = @(td, yd) Desorption(td, yd, m_H2_out(i), Tw_in_des(i), DTml_des(i), C_des, E_des, R, P_des, m_s, M_H2, SC, M_MH, Cp_H2, Cp_s, A, U, DH_des, DS_des, sl, P0);
    
    % Differential equation solver
    [td, yd] = ode23s(ode_fun_d, tspan, y0_d);

    % Store solutions for current timestep
    m_H2_des(i+1) = yd(end, 1);
    m_MH_des(i+1) = yd(end, 2);
    T_des(i+1)    = yd(end, 3);

    % Calculations of the other parameters
    Peq_des(i) = (exp(-(DH_des/(R*T_des(i)))+(DS_des/R)+sl*((m_MH_des(i)/m_s)-0.5)))*P0;           % Desorption equilibrium pressure           [Pa]
    r_des(i)   = C_des*exp(-E_des/(R*T_des(i)))*((P_des-Peq_des(i))/Peq_des(i))*(m_MH_des(i)/m_s); % Reaction rate                             [1/s]

    % Calculate outlet water temperature
    DTb_des(i+1)  = Tw_in_des(i) - T_des(i);                                                       % Calculation of the new DTb                [K]
    DTa_des(i+1)  = Tw_out_des(i) - T_des(i);                                                      % Calculation of the new DTa                [K]
    DTml_des(i+1) = (DTb_des(i+1) - DTa_des(i+1)) / (log(DTb_des(i+1) / DTa_des(i+1)));            % Calculation of the new DTml               [K]
    Tw_out_des(i+1) = Tw_in_des(i) - (DTml_des(i+1)*U*A)/(m_w_des(i)*Cp_w*1000);                   % Calculation of the new Tw_out             [K]
    Q_des(i+1)        = (m_w_des(i)*Cp_w*(Tw_in_des(i)-Tw_out_des(i+1)));                          % Calculation of the heat (water --> MH)    [kW]

    res_des       = DTml_des(i+1)*U*A/1000 - (Tw_in_des(i)-Tw_out_des(i+1))*m_w_des(i)*Cp_w;       % Energy balance equation to be verified
    end
end

Qdes = Q_des'; % Transposition of the row vector in a column vector 

%% Plot result
% Temperatures
figure(1)

    %subplot(2,1,[1 1]);
    plot(t_abs, T_abs,'k',t_des, T_des,'r','LineWidth',1.5);
    ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    ylim([260 360]);
    title('Metal Hydride Temperature "T_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    xlabel('Time [s]','Fontsize',25, 'FontName','Times New Roman','FontWeight','bold');
    xlim([t_abs(1)-100 t_des(end)+100]);
    set(gca,'FontSize',25,'FontAngle','Normal','FontName','Times New Roman');
    legend('Absorption', 'Desorption','location','northwest','FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    hold on
    grid on

    % subplot(2,2,3);
    % riga1 = 'Absorption';
    % riga2 = 'Outlet water temperature "T_{w,out}"';
    % titolo = sprintf('%s\n%s', riga1, riga2);
    % plot(t_abs, Tw_out_a,'k');
    % ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % ylim([280 295]);
    % xlim([-100 2000]);
    % xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % title(titolo, 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % hold on
    % grid on
    % 
    % subplot(2,2,4)
    % riga1 = 'Desorption';
    % riga2 = 'Outlet water temperature "T_{w,out}"';
    % titolo = sprintf('%s\n%s', riga1, riga2);
    % plot(t_des, Tw_out_d,'r');
    % ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % ylim([330 355]);
    % title(titolo, 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % xlim([1700 3700]);
    % xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    % hold on
    % grid on
 

figure (2)

% Hydride metal mass
    subplot(2, 2, 1);
    plot(t_abs,m_MH_abs,'k',t_des, m_MH_des,'r','LineWidth',1.5);
    ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    title('Metal Hydride mass "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    xlim([t_abs(1)-100 t_des(end)+100]);
    set(gca,'FontSize',25,'FontAngle','Normal','FontName','Times New Roman');
    hold on
    grid on

% Equilibrium pressure
    subplot(2, 2, 2);
    plot(t_abs,Peq_abs/1E+05,'k',t_des, Peq_des/1E+05,'r','LineWidth',1.5);
    ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    title('Equilibrium pressure "P_{eq}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    xlim([t_abs(1)-100 t_des(end)+100]);
    set(gca,'FontSize',25,'FontAngle','Normal','FontName','Times New Roman');
    hold on
    grid on

% Reaction rate
    subplot(2, 2, 3);
    plot(t_abs, r_abs,'k',t_des, r_des,'r','LineWidth',1.5);
    ylabel('[1/s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    title('Reaction rate "r"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    xlim([t_abs(1)-100 t_des(end)+100]);
    set(gca,'FontSize',25,'FontAngle','Normal','FontName','Times New Roman');
    hold on
    grid on

% Hydrogen mass
    subplot(2, 2, 4);
    plot(t_abs, m_H2_abs,'k',t_des, m_H2_des,'r','LineWidth',1.5);
    ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',25);
    title('Hydrogen mass "m_{H_2}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',18);
    xlim([t_abs(1)-100 t_des(end)+100]);
    set(gca,'FontSize',25,'FontAngle','Normal','FontName','Times New Roman');
    hold on
    grid on

toc
