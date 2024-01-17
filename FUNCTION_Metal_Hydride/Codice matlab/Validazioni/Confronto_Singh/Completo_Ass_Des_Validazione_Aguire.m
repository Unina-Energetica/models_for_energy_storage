clear all; clc;

%% Equations to be used
% Calculation of m_s based on the quantity of m_H2 to be absorbed
% m_s = (m_H2/wt_max) * 100  [kg]  con wt_max = 1.39
%
% Calculation of hydrogen flowrate m_H2_in based on the quantity absorbed
% m_H2_in = (2 * m_H2) / 3600  [kg/s]

%% Parameters absorption and desorption
rho_MH = 8310;        % Density of LaNi5                   [kg/m3]
Vg     = 0.0172;      % Gas Volume                         [m3]
R      = 8.314;       % Gas constant                       [J/mol K]
M_H2   = 2.016/1000;  % Molar mass H2                      [kg/mol]
M_MH   = 0.432;       % Molar mass   MH                    [kg/mol]
Cp_H2  = 14890;       % Specific heat H2                   [J/kg K]
Cp_s   = 385;         % Specific heat solid                [J/kg K]
A      = 5;         % Area "heat interchange"            [m2]
U      = 325;         % Overall heat transfer coefficient  [W/m2 K]
sl     = 0.13;        % Slope coefficient                  [-]
P0     = 1E+05;       % Reference pressure                 [Pa]
m_s    = 143*2.9;         % Solid mass                         [kg]
SC     = 2.87;           % Stoichiometric coefficient         [-]
Cp_w   = 4.185;       % Specific heat water                [kJ/kg K] 

%% Simulation time for absorption and desorption
t_abs  = (0:1:3600)';              % [s]
t_des  = (t_abs(end)+1:1:2*t_abs(end)+1)';   % [s]

%% Input parameters
m_H2_in  = ones(numel(t_abs),1)*(2*10/3600);   % Input hydrogen flowrate               [kg/s]
% m_H2_in= xlsread('m_H2_in_TRNSYS.xlsx','Foglio1','A1:A1801');
m_H2_out = ones(numel(t_des),1)*(2*4/3600);   % Output hydrogen flowrate              [kg/s]
V_Twa    = ones(numel(t_abs),1) * 300;      % Inlet cooling water temperature       [k]
%V_Twa   = rand(numel(t_abs),1).*10+Twa;  
V_Twd    = ones(numel(t_des),1) * 353;      % Inlet heating water temperature fixed [K]
m_w_abs  = ones(numel(t_abs),1) * 0.947;   % Project water flowrate for absor      [kg/s]
m_w_des  = ones(numel(t_des),1) * 1.1947;   % Project water for desorption          [kg/s]

%% Absorption
%Absorption parameters
Ca      = 59.187;   % Kinetic constant                   [1/s]
Ea      = 21170;  % Absorption activation energy       [J/mol]
Pa      = 15E+05;  % Absorption pressure                [Pa]
Tin_a   = 300;    % H2 input temperature               [K]
DHa     = -30478; % Absorption enthalpy                [J/mol]
DSa     = -108;   % Absorption entropy                 [J/mol K]
eps_abs = 10E-10; % Control value for energy balance   [-]

% Variable inizialization 
m_H2_a   = zeros(size(t_abs));
m_MH_a   = zeros(size(t_abs));
T_a      = zeros(size(t_abs));
Tw_out_a = zeros(size(t_abs));
DTb_abs  = zeros(size(t_abs));
DTa_abs  = zeros(size(t_abs));
DTml_abs = zeros(size(t_abs));
Peq_a    = zeros(size(t_abs));
r_a      = zeros(size(t_abs));

% Initial conditions for the first timestep
m_H2_a(1)   = 0;      % Absorbed hydrogen mass           [kg]
m_MH_a(1)   = 0;      % Metal hydride mass               [kg]
T_a(1)      = V_Twa(1)+6;    % System's temperature             [K]
Tw_out_a(1) = V_Twa(1)+5; % Initial outlet water temperature [K] 
DTb_abs(1)  = T_a(1) - V_Twa(1);
DTa_abs(1)  = T_a(1) - Tw_out_a(1);
DTml_abs(1) = (DTb_abs(1)-DTa_abs(1))/(log(DTb_abs(1)/DTa_abs(1)));

for i = 1:numel(t_abs)-1    % We use "for" in order to have 1801 iterations

    res_abs = 1000;     % Initial value for the residual

    while res_abs > eps_abs % We repeat the iteration if the energy balance is not verified

    if i == 1           % We use "if" to assign the initial condition based on the "i" index
        y0_a        = [m_H2_a(1); m_MH_a(1); T_a(1)];
        DTb_abs(i)  = DTb_abs(1);
        DTa_abs(i)  = DTa_abs(1);
        DTml_abs(i) = (DTb_abs(1)-DTa_abs(1))/(log(DTb_abs(1)/DTa_abs(1)));
    else
        y0_a        = [m_H2_a(i); m_MH_a(i); T_a(i)];
    end

    % Define tspan for ODE resolution for current timestep
    tspan = [t_abs(i) t_abs(i+1)];

    % Function handle
    ode_fun_a = @(ta, ya) Assorbimento_for(ta, ya, m_H2_in(i), V_Twa(i), DTml_abs(i), m_s, A, U, Pa, Ca, Ea, R, M_H2, SC, M_MH, Cp_H2, Cp_s, Tin_a, DHa, DSa, sl, P0);

    % Differential equation solver
    [ta, ya] = ode23s(ode_fun_a, tspan, y0_a);

    % Store solutions for current timestep 
    m_H2_a(i+1) = ya(end, 1);
    m_MH_a(i+1) = ya(end, 2);
    T_a(i+1)    = ya(end, 3);

    % Calculation of other parameters
    Peq_a(i+1) = exp((DHa/(R*T_a(i)))-DSa/R+sl*((m_MH_a(i)/m_s)-0.5))*P0;   % Absorption equilibrium pressure [Pa]
    r_a(i+1)   = Ca*exp(-Ea/(R*T_a(i)))*log(Pa/Peq_a(i))*(1-(m_MH_a(i)/m_s));  % Reaction rate                   [1/s]
    wt(i+1)    = (m_MH_a(i)/m_s)*(M_H2*SC/M_MH)*100;

    % Calculate outlet water temperature
    DTb_abs(i+1)    = T_a(i) - V_Twa(i);                                                 % Calculation of the new DTb
    DTa_abs(i+1)    = T_a(i) - Tw_out_a(i);                                              % Calculation of the new DTa
    DTml_abs(i+1)   = (DTb_abs(i+1) - DTa_abs(i+1)) / log(DTb_abs(i+1) / DTa_abs(i+1));  % Calculation of the new DTml
    Tw_out_a(i+1)   = V_Twa(i) + (DTml_abs(i+1) * U * A) / (m_w_abs(i) * Cp_w * 1000);   % Calculation of the new Tw_out
    Q(i+1)          = (m_w_abs(i)*Cp_w*(Tw_out_a(i+1)-V_Twa(i)));                        % Calculation of the heat (MH --> water)

    res_abs         = DTml_abs(i+1)*U*A/1000 - (Tw_out_a(i+1)-V_Twa(i))*m_w_abs(i)*Cp_w; % Energy balance equation to be verified 
    end
end

Q = Q'; % Transposition of the row vector in a column vector 
wt = wt';

%% Desorption
% Desorption parameters

Cd      = 9.57;      % Kinetic constant              [1/s]
Ed      = 16450;    % Desorption activation energy  [J/mol]
Pd      = 6E+05;    % Desorption pressure           [Pa]
DHd     = 30800;    % Desorption enthalpy           [J/mol]
DSd     = 108;      % Desorption entropy            [J/mol K]
eps_des = 10E-10;   % Control value for energy balance

% Variable inizialization
m_H2_d   = zeros(size(t_des));
m_MH_d   = zeros(size(t_des));
T_d      = zeros(size(t_des));
Tw_out_d = zeros(size(t_des));
DTb_des  = zeros(size(t_des));
DTa_des  = zeros(size(t_des));
DTml_des = zeros(size(t_des));
Peq_d    = zeros(size(t_des));
r_d      = zeros(size(t_des));

% Initial conditions for the first timestep
m_H2_d(1)   = 0;
m_MH_d(1)   = m_MH_a(end);
T_d(1)      = T_a(end);
Tw_out_d(1) = V_Twd(1)-10;
DTb_des(1)  = V_Twd(1) - T_d(1);
DTa_des(1)  = Tw_out_d(1) - T_d(1);
DTml_des(1) = (DTb_des(1)-DTa_des(1))/(log(DTb_des(1)/DTa_des(1)));
Peq_d(1)    = Peq_a(end);


for i = 1:numel(t_des)-1

    res_des = 1000;

    while res_des > eps_des

        if i == 1
            y0_d        = [m_H2_d(1); m_MH_d(1); T_d(1)];
            DTb_des(i)  = DTb_des(1);
            DTa_des(i)  = DTa_des(1);
            DTml_des(i) = DTml_des(1);
        else
            y0_d = [m_H2_d(i); m_MH_d(i); T_d(i)];
        end

    % Define tspan for ODE resolution for current timestep
    tspan = [t_des(i) t_des(i+1)];

    % Function handle
    ode_fun_d = @(td, yd) Desorbimento_for(td, yd, m_H2_out(i), V_Twd(i), DTml_des(i), Cd, Ed, R, Pd, m_s, M_H2, SC, M_MH, Cp_H2, Cp_s, A, U, DHd, DSd, sl, P0);
    
    % Differential equation solver
    [td, yd] = ode23s(ode_fun_d, tspan, y0_d);

    % Store solutions for current timestep
    m_H2_d(i+1) = yd(end, 1);
    m_MH_d(i+1) = yd(end, 2);
    T_d(i+1)    = yd(end, 3);

    % Calculations of the other parameters
    Peq_d(i) = (exp(-(DHd/(R*T_d(i)))+(DSd/R)+sl*((m_MH_d(i)/m_s)-0.5)))*P0;    % Desorption equilibrium pressure [Pa]
    r_d(i)   = Cd*exp(-Ed/(R*T_d(i)))*((Pd-Peq_d(i))/Peq_d(i))*(m_MH_d(i)/m_s); % Reaction rate

    % Calculate outlet water temperature
    DTb_des(i+1)  = V_Twd(i) - T_d(i);                                                  % Calculation of the new DTb
    DTa_des(i+1)  = Tw_out_d(i) - T_d(i);                                               % Calculation of the new DTa
    DTml_des(i+1) = (DTb_des(i+1) - DTa_des(i+1)) / (log(DTb_des(i+1) / DTa_des(i+1))); % Calculation of the new DTml
    Tw_out_d(i+1) = V_Twd(i) - (DTml_des(i+1)*U*A)/(m_w_des(i)*Cp_w*1000);              % Calculation of the new Tw_out
    Q(i+1)        = (m_w_des(i)*Cp_w*(V_Twd(i)-Tw_out_d(i+1)));                         % Calculation of the new heat power (water --> MH)

    res_des       = DTml_des(i+1)*U*A/1000 - (V_Twd(i)-Tw_out_d(i+1))*m_w_des(i)*Cp_w;  % Energy balance equation to be verified
    end
end


% %% Plot result
% 
% % Temperatures
figure(1)

    subplot(2,1,[1 1]);
    plot(t_abs, T_a,'k',t_des, T_d,'r');
    ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    ylim([260 360]);
    title('Metal Hydride Temperature "T_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([t_abs(1)-100 t_des(end)+100]);
    legend('Absorption', 'Desorption','location','northwest','FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

    subplot(2,2,3);
    riga1 = 'Absorption';
    riga2 = 'Outlet water temperature "T_{w,out}"';
    titolo = sprintf('%s\n%s', riga1, riga2);
    plot(t_abs, Tw_out_a,'k');
    ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    ylim([280 295]);
    xlim([-100 2000]);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title(titolo, 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on

    subplot(2,2,4)
    riga1 = 'Desorption';
    riga2 = 'Outlet water temperature "T_{w,out}"';
    titolo = sprintf('%s\n%s', riga1, riga2);
    plot(t_des, Tw_out_d,'r');
    ylabel('[K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    ylim([330 355]);
    title(titolo, 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([1700 3700]);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on
    grid on


figure (2)

% Hydride metal mass
    subplot(2, 2, 1);
    plot(t_abs,m_MH_a,'k',t_des, m_MH_d,'r');
    ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Metal Hydride mass "m_{MH}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([t_abs(1)-100 t_des(end)+100]);
    hold on
    grid on

% Equilibrium pressure
    subplot(2, 2, 2);
    plot(t_abs,Peq_a/1E+05,'k',t_des, Peq_d/1E+05,'r');
    ylabel('[bar]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Equilibrium pressure "P_{eq}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([t_abs(1)-100 t_des(end)+100]);
    hold on
    grid on

% Reaction rate
    subplot(2, 2, 3);
    plot(t_abs, r_a,'k',t_des, r_d,'r');
    ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Reaction rate "r"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([t_abs(1)-100 t_des(end)+100]);
    hold on
    grid on

% Hydrogen mass
    subplot(2, 2, 4);
    plot(t_abs, m_H2_a,'k',t_des, m_H2_d,'r');
    ylabel('[kg]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Hydrogen mass "m_{H_2}"', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    xlim([t_abs(1)-100 t_des(end)+100]);
    hold on
    grid on


   %% Risultati paper
   wt_paper_15bar = xlsread('wttt.xlsx','Foglio1','B1:B26');
   t_15bar_wt_pap = xlsread('wttt.xlsx','Foglio1','A1:A26');
   wt_paper_5bar  = xlsread('wttt.xlsx','Foglio1','D1:D32');
   t_5bar_wt_pap  = xlsread('wttt.xlsx','Foglio1','C1:C32');

   wt_valid_15bar = interp1(t_abs, wt, t_15bar_wt_pap, 'linear');
   wt_valid_5bar  = interp1(t_abs, wt, t_5bar_wt_pap, 'linear');

   if Pa == 15E+05
        figure(4)
        subplot(1,2,1);
         yyaxis left
         plot(t_15bar_wt_pap, wt_paper_15bar, 'ks');
         ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         hold on
         yyaxis right 
         plot(t_15bar_wt_pap, wt_valid_15bar);
         ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         hold on
         xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         title('Comparison model vs experimental', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            legend('Riferimento', 'Modello UNINA');

   elseif Pa == 5E+05
        figure(4)
        subplot(1,2,1);
         yyaxis left
         plot(t_5bar_wt_pap, wt_paper_5bar, 'ks');
         ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         hold on
         yyaxis right 
         plot(t_5bar_wt_pap, wt_valid_5bar);
         ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         hold on
         xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
         title('Comparison model vs experimental', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            legend('Riferimento', 'Modello UNINA');

   end


   t_paper_15bar = xlsread('Temp.xlsx','Foglio1','A1:A44');
   T_paper_15bar = xlsread('Temp.xlsx','Foglio1','B1:B44');
   t_paper_5bar  = xlsread('Temp.xlsx','Foglio1','C1:C31');
   T_paper_5bar  = xlsread('Temp.xlsx','Foglio1','D1:D31');

   T_valid_15bar = interp1(t_abs, T_a, t_paper_15bar, 'linear');
   T_valid_5bar  = interp1(t_abs, T_a, t_paper_5bar, 'linear');

      if Pa == 15E+05
        figure(4)
        subplot(1,2,2);
            yyaxis left
            plot(t_paper_15bar, T_paper_15bar, 'ks');
            ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            ylim([290 360]);
            hold on
            yyaxis right 
            plot(t_paper_15bar, T_valid_15bar);
            ylim([290 360]);
            ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            hold on
            xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            title('Comparison model vs experimental', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            legend('Riferimento', 'Modello UNINA');
      elseif Pa == 5E+05
          figure(4)
          subplot(1,2,2);
            yyaxis left
            plot(t_paper_5bar, T_paper_5bar, 'ks');
            ylim([290 360]);
            ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            hold on
            yyaxis right 
            plot(t_paper_5bar, T_valid_5bar);
            ylim([290 360]);
            ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            hold on
            xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            title('Comparison model vs experimental', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
            legend('Riferimento', 'Modello UNINA');
      end

if Pa == 15E+05
    figure(5)
    yyaxis left
    plot(t_15bar_wt_pap, wt_paper_15bar, 'ko', t_15bar_wt_pap, wt_valid_15bar);
    ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    ylim([0 1.5]);
        legend('Riferimento', 'Modello UNINA');
    yyaxis right
    plot(t_paper_15bar, T_paper_15bar, 'ks', t_paper_15bar, T_valid_15bar);
    ylim([290 360]);
    ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on 
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Storage capacity and temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legend('Riferimento', 'Modello UNINA');

elseif Pa == 5E+05
        figure(5)
    yyaxis left
    plot(t_5bar_wt_pap, wt_paper_5bar, 'ko', t_5bar_wt_pap, wt_valid_5bar);
    ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    ylim([0 1.5]);
        legend('Riferimento', 'Modello UNINA');
    yyaxis right
    plot(t_paper_5bar, T_paper_5bar, 'ks', t_paper_5bar, T_valid_5bar);
    ylim([290 360]);
    ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    hold on 
    xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    title('Storage capacity and temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
    legend('Riferimento', 'Modello UNINA');
end


% figure(6)    
%     yyaxis left
%     plot(t_abs, wt);
%     ylabel('Storage capacity [%]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     ylim([0 1.5]);
%     yyaxis right
%     plot(t_abs, T_a-273.15);
%     ylim([10 70]);
%     ylabel('Temperature [K]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     hold on 
%     xlabel('Time [s]', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);
%     title('Storage capacity and temperature', 'FontName','Times New Roman','FontWeight','bold', 'FontSize',13);

% Error calculation
diff_perc_T = ((T_paper_15bar - T_valid_15bar) ./ T_paper_15bar) * 100;
% diff_perc_wt = ((wt_paper_15bar - wt_valid_15bar) ./ wt_paper_15bar) * 100;

figure
plot(t_paper_15bar, diff_perc_T);
ylim([-10 10]);

%figure
% plot(t_15bar_wt_pap, diff_perc_wt);

