function [T_gas_out1,Ts_out,Qheat_sun,Eps_sun,T_gas_out2,To_out,Qheat_oil,Eps_oil]=heater3V2(T_gas_in,Ts_in,To_in,m_gas,m_heater_sun,m_heater_oil)
%% HEATER 3 %% 

% clear all
% T_gas_in = 30;
% Ts_in = 120;
% To_in = 210;
% m_gas = 3;
% m_heater_sun = 5.5;
% m_heater_oil = 5.5;

%% Condizioni di progetto aria-olio (solare)

T_air_in1 = 30; % °C
m_air = 3.4; % kg/s

T_sun_in = 110;
% m_sun = 6;
T_sun_out = 90;

DTs = T_sun_in - T_sun_out;

delta_Ta = 10;

cp_air = 1.0204; % kJ/ kg K
cp_sun = 2.195; % kJ/ kg K

T_air_out1 = T_sun_in - delta_Ta;

Q_rif = m_air*cp_air*(T_air_out1 - T_air_in1);

% T_sun_out = T_sun_in - Q_rif/(cp_sun*m_sun);
m_sun = Q_rif/(DTs*cp_sun);

Cmin_ref = min(cp_sun*m_sun,cp_air*m_air);
Cmax_ref = max(cp_sun*m_sun,cp_air*m_air);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in1 - T_sun_out);
DTb = -(T_air_out1 - T_sun_in);
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %W/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_rif = UA_rif/Cmin_ref;
% return

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_sun_new = cp_sun*m_heater_sun; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_sun_new);
Cmax_new = max(C_air_new,C_sun_new);

U_rap = (m_gas/m_air)^0.8;

NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = -(T_gas_in - Ts_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out1 = T_gas_in + Q_int/(cp_air*m_gas); % kW (K/kW)
Ts_out = Ts_in - Q_int/(cp_sun*m_heater_sun);

Qheat_sun  = Q_int;
Eps_sun    = epsilon_int;

%% Condizioni di progetto aria-olio 

T_air_in2 = T_air_out1; % °C
% m_air = 3.2; % kg/s

T_oil_in = 250;
% m_oil = 6.4;
T_oil_out = 210;

DTo = T_oil_in - T_oil_out;

% delta_Ta = 10;

% cp_air = 1.4404; % kJ/ kg K
cp_oil = 2.195; % kJ/ kg K

T_air_out2 = 240;

Q_rif = m_air*cp_air*(T_air_out2 - T_air_in2);

% T_oil_out = T_oil_in - Q_rif/(cp_oil*m_oil);
m_oil = Q_rif/(DTo*cp_oil);

Cmin_ref = min(cp_oil*m_oil,cp_air*m_air);
Cmax_ref = max(cp_oil*m_oil,cp_air*m_air);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in2 - T_oil_out);
DTb = -(T_air_out2 - T_oil_in);
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %W/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_rif = UA_rif/Cmin_ref;
% return

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_oil_new = cp_oil*m_heater_oil; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_oil_new);
Cmax_new = max(C_air_new,C_oil_new);

U_rap = (m_gas/m_air)^0.8;

NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = -(T_gas_out1 - To_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out2 = T_gas_out1 + Q_int/(cp_air*m_gas);  % kW (K/kW)
To_out = To_in - Q_int/(cp_oil*m_heater_oil);

Qheat_oil  = Q_int;
Eps_oil    = epsilon_int;

