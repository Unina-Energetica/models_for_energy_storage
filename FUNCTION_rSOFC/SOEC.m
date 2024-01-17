% SOEC model
clear; close all; clc

% Parameters

R         = 8.314;    % Gas constant (m3*Pa/mol*K) or (J/mol*K)
F         = 96485;    % Faraday constant (C/mol)

p_anode   = 1;        % anode pressure   (bar)
p_cathode = 1;        % cathode pressure (bar)
Acell     = 87.7;     % cell active area (cm2)
ncell     = 8;        % number of cells 
nstack    = 2500;     % number of stacks 
Uf        = 0.8;      % fuel utilization factor

% upper and lower limit are selected to allow operation of the cell in the 
% optimal range (300 - 1100 mA/cm2)

% Input variables

P_SOEC=30*nstack:nstack/5:1140*nstack;

V_ideal   = 10;                                          % reference voltage of the stack (V) 
I_SOEC    = P_SOEC/V_ideal;                              % operating current range (A)

for i=1:length(I_SOEC) 
    
if I_SOEC(i) == 0
 
    W_H2(i)     = 0;                                     % H2 volumetric flow rate (m3/h)
    T           = 0;                                     % cell temperature (K)
    eta(i)      = 0;                                     % cell efficiency 
    Vstack(i)   = 0;                                     % stack voltage (V)
    Pel(i)      = 0;                                     % input electric power (W)
    Qheat(i)    = 0;                                     % input heat (W)

else

    Istack(i)   = I_SOEC(i)./nstack;                     % stack current (A)
    Icell(i)    = Istack(i);                             % cell current (A)
    jcell(i)    = (Icell(i)./Acell).*1000;               % cell density current (mA/cm2)
    z_dot(i)    = (Icell(i))./(2*F);                     % H2O molar flow rate tacking part to the reaction (mol/s) 
                                             
    nH2O_in(i)  = z_dot(i)./Uf;                          % H2O input molar flow rate (mol/s)
    xH2O_in     = 1;                                     % H2O input molar fraction (mol/s)
    nH2_out(i)  = z_dot(i);                              % H2 output molar flow rate (mol/s)
    nO2_out(i)  = 0.5.*z_dot(i);                         % O2 output molar flow rate (mol/s)
    nH2O_out(i) = nH2O_in(i)-z_dot(i);                   % H2O output molar flow rate (mol/s)
    xO2_out     = 1;                                     % O2 output molar fraction
    xH2O_out(i) = nH2O_out(i)./(nH2_out(i)+nH2O_out(i)); % H2O output molar fraction
    xH2_out(i)  = nH2_out(i)./(nH2_out(i)+nH2O_out(i));  % H2 output molar fraction

    % Enthalpy
    T0          = 273.15;                                % standard temperature condition (K)
    Tin         = 1073;                                  % inlet fuel temperature (K)
    A_H2O = 3.47;    A_H2 = 3.3249;   A_O2 = 3.639;
    B_H2O = 0.00145; B_H2 = 0.000422; B_O2 = 0.000506;
    D_H2O = 1210;    D_H2 = 8300;     D_O2 = -22700;
    
    % IN
    hH2O_in     = R.*(A_H2O.*(Tin-T0)+B_H2O./2.*((Tin.^2)-(T0.^2))-D_H2O.*(1./Tin-1./T0));  % (J/mol)
    % OUT
    hH2_out     = @(T) R.*(A_H2.*(T-T0)+B_H2./2.*((T.^2)-(T0.^2))-D_H2.*(1./T-1./T0));      % (J/mol)
    hO2_out     = @(T) R.*(A_O2.*(T-T0)+B_O2./2.*((T.^2)-(T0.^2))-D_O2.*(1./T-1./T0));      % (J/mol)
    hH2O_out    = @(T) R.*(A_H2O.*(T-T0)+B_H2O./2.*((T.^2)-(T0.^2))-D_H2O.*(1./T-1./T0));   % (J/mol)
    
    % Gibbs free energy
    T0          = 298.15;                                % standard temperature condition (K)
    gf          = -228572;                               % standard formation Gibbs free energy (J/mol)
    hf          = -241818;                               % standard formation enthalpy (J/mol)

    A  = A_H2O-A_H2-0.5*A_O2; B  = B_H2O-B_H2-0.5*B_O2; C  = 0; D  = D_H2O-D_H2-0.5*D_O2;

    g           = @(T) (gf-hf)./(R.*T0)+hf./(R.*T)+(T.*A+1/2*T.^2.*B+1/3.*T.^3.*C-1./T.*D-T0.*A-1/2.*T0.^2.*B-1/3.*T0.^3.*C+1./T0.*D)./T-(1/2*C.*T.^2+B.*T-1/2*D./T.^2+A.*log(T)-1/2*T0.^2*C-T0.*B+1/2.*D./T0.^2-A.*log(T0));
    gcell       = @(T) g(T).*R.*T;   % (J/mol)
    % theoretical efficiency
    etamax      = @(T) gcell(T)./hf;                    
    % open circuit voltage at standard condition
    E0          = @(T) -gcell(T)./(2*F);
    % open circuit voltage
    Ecell       = @(T) E0(T)+R.*T.*log((xH2_out(i).*p_anode).*(xO2_out.*p_cathode).^0.5./(xH2O_in.*p_anode))./(2*F); 
  
    % Activation losses
    g_a         = 1344000;                                         % (A/cm2)
    g_c         = 205100;                                          % (A/cm2)
    ea_a        = 100000;                                          % (J/mol)
    ea_c        = 120000;                                          % (J/mol)
    i0a         = @(T) g_a.*exp(-ea_a./(R.*T))*1000;               % (mA/cm2)
    i0c         = @(T) g_c.*exp(-ea_c./(R.*T))*1000;               % (mA/cm2)
    Vacta       = @(T) R.*T./F.*asinh(jcell(i)./(2.*i0a(T)));                                        
    Vactc       = @(T) R.*T./(2*F).*asinh(jcell(i)./(2.*i0c(T)));                                    
    Vact        = @(T) Vacta(T) + Vactc(T);                        % overvoltage due to activation losses     
   
    % Concentration losses
    il          = 1340;                                            % upper limit current density (mA/cm2)
    Vconc       = @(T) -R*T/(2*F)*log((1-jcell(i)./il));           % overvoltage due to concentration losses  
    
    % Ohmic losses
    sigma0      = 333.3;                                           % Pre-exponential factor (1/(ohm*cm))
    ea_el       = 85634;                                           % Electrolyte activation energy (J/mol)
    del         = 0.00125;                                         % Electrolyte thickness (cm)  
    sigma_el    = @(T) sigma0.*exp(-ea_el./(R.*T));                % (1/(ohm*cm))
    rel         = @(T) del./sigma_el(T);                           % Electrolyte ohmic resistance (ohm*cm2)
    rstack      = 0.057;                                           % Stack ohmic resistance (ohm*cm2)
    rtot        = @(T) (rel(T) + rstack);                          % Total ohmic resistance (ohm*cm2)
    Vohm        = @(T) rtot(T).*(jcell(i)./1000);                  % overvoltage due to ohmic losses (V)

    % Total cell voltage
    Vcell       = @(T) Ecell(T) + Vact(T) + Vconc(T) + Vohm(T); 
    
    % Balance
    f1          = (nH2O_in(i).*hH2O_in);                           % W
    f2          = @(T) -nH2_out(i).*hH2_out(T);                    % W
    f3          = @(T) -nH2O_out(i).*hH2O_out(T);                  % W
    f4          = @(T) -nO2_out(i).*hO2_out(T);                    % W
    f5          = @(T) Vcell(T).*(jcell(i)./1000).*Acell;          % W                               
    f6          = (z_dot(i).*hf);                                  % W                             
    energy      = @(T) f1 + f2(T) + f3(T) + f4(T) + f5(T) + f6;    % W
    
    % fzero
    T = fzero(energy,1000);
    
        % Cell Energy (W)
        Energy(i)    = energy(T);              
        % Cell current (A)
        Icell(i)     = jcell(i)./1000.*Acell;    
        % Cell voltage (V)
        V_cell(i)    = Vcell(T);               
        E_cell(i)    = Ecell(T);               
        V_act(i)     = Vact(T);                 
        V_conc(i)    = Vconc(T);              
        V_ohm(i)     = Vohm(T);               
        % Cell power (W) 
        pcell(i)     = Icell(i)*Vcell(T);       
        % Stack power (W)
        Pstack(i)    = pcell(i)*ncell;         
        % Stack voltage (V)
        Vstack(i)    = Vcell(T).*ncell;      
        % SOEC input power (W)
        Pel(i)       = Pstack(i)*nstack;       
        % Q heat H2O (W)
        Tin_cella    = 1073; 
        T0_cella     = 298.15;
        hH2O         = R.*(A_H2O.*(Tin_cella-T0_cella)+B_H2O./2.*((Tin_cella.^2)-(T0_cella.^2))-D_H2O.*(1./Tin_cella-1./T0_cella));  
        Qheat_H2O(i) = nH2O_in(i).*hH2O;          
        Qheat(i)     = Qheat_H2O(i).*ncell*nstack;

        % SOEC efficiency
        LHV_H2       = 120;                                                        % Lower Heating Value (MJ/kg)
        PM_H2        = 2.016;                                                      % Molar weight (g/mol)
        mH2_out(i)   = nH2_out(i).*PM_H2;                                          % Mass flow rate (g/s)
        eta(i)       = (mH2_out(i).*LHV_H2*1000/(pcell(i)+Qheat_H2O(i)))*100;      % Efficiency 
        W_H2(i)      = (nH2_out(i).*ncell*nstack*3600)*(R/100000)*T/p_anode;       % Volumetric flow rate (m3/h)
  
end
    ndot_H2_tot(i) = nH2_out(i).*ncell*nstack*3600;        % Total hourly molar flow rate (mol/h)
    mdot_H2_tot(i) = mH2_out(i).*ncell*nstack*3600/1000;   % Total hourly mass flow rate (kg/h)

RES.m_H2(i)  = mdot_H2_tot(i);
RES.n_H2(i)  = ndot_H2_tot(i);
RES.Pel(i)   = Pel(i)/1000;
RES.eta(i)   = eta(i);
RES.V(i)     = Vstack(i);
RES.T(i)     = T;

end


RES.m_H2(end)
RES.n_H2(end)
RES.Pel(end)
RES.eta(end)
RES.T(end)

jcell_start  = jcell(1);
jcell_end    = jcell(end);
I_SOEC_start = I_SOEC(1);
I_SOEC_end   = I_SOEC(end);

figure(1)
plot(jcell,RES.T)
xlabel('density current [mA/cm2]')
ylabel('cell temperature [K]')
xlim([jcell_start jcell_end])
grid on
legend('T')

figure(2)
plot(jcell,V_cell,'r')
title('Polarization curve')
xlabel('density current [mA/cm2]')
ylabel('cell voltage [V]')
xlim([jcell_start jcell_end])
grid on
legend('V')

figure(3)
% plot(jcell,RES.m_H2*1000)
% hold on
plot(jcell,RES.Pel)
hold on
xlabel('density current [mA/cm2]')
ylabel('SOEC power [kW]')
yyaxis('right')
ylabel('\eta_{SOEC}')
plot(jcell,RES.eta)
xlim([jcell_start jcell_end])
grid on
% legend('m_{H2} [mg/h]','Pel [kW]','\eta_{SOEC}')
legend('Pel [kW]','\eta_{SOEC}')

figure(4)
plot(jcell,V_cell)
hold on
plot(jcell,E_cell)
hold on
plot(jcell,V_conc)
hold on
plot(jcell,V_ohm)
hold on 
plot(jcell,V_act)
xlabel('density current [mA/cm^{2}]')
ylabel('cell voltage [V]')
xlim([jcell_start jcell_end])
grid on
legend('Vcell','Ecell','Vconc','Vohm','Vact')

figure(5)
plot(jcell,Vstack)
xlabel('density current [mA/cm^{2}]')
ylabel('SOEC voltage [V]')
xlim([jcell_start jcell_end])
grid on
