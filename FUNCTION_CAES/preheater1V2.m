function [T_gas_out1,Tamb_out,Qheat_amb,Eps_amb,T_gas_out2,Tw_out,Qheat_water,Eps_water,T_gas_out3,Ts_out,Qheat_sun,Eps_sun,T_gas_out4,To_out,Qheat_oil,Eps_oil]=preheater1V2(T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,m_gas,m_heater_amb,m_heater_water,m_heater_sun,m_heater_oil)
%% PREHEATER 1 %% 

% clear all
% T_gas_in = 200-273.15;
% Tamb_in = 18;
% Tw_in = 31;
% Ts_in = 120;
% To_in = 215;
% m_gas = 2;
% m_heater_amb = 30;
% m_heater_water = 8;
% m_heater_sun = 9;
% m_heater_oil = 13.4;

%% Condizioni di progetto aria-aria 

T_air_in1 = 180-273.15; % °C
m_air = 3.4; % kg/s

T_amb_in = 20;
m_amb = 30;

delta_Ta = 10;

cp_air = 1.4404; % kJ/ kg K
cp_amb = 1.006; % kJ/ kg K

Cmin_ref = min(cp_amb*m_amb,cp_air*m_air);
Cmax_ref = max(cp_amb*m_amb,cp_air*m_air);

T_air_out1 = T_amb_in - delta_Ta;

Q_rif = m_air*cp_air*(T_air_out1 - T_air_in1);

T_amb_out = T_amb_in - Q_rif/(cp_amb*m_amb);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in1 - T_amb_out);
DTb = -(T_air_out1 - T_amb_in);
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %W/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_int = UA_rif/Cmin_ref;
% return

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_amb_new = cp_amb*m_heater_amb; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_amb_new);
Cmax_new = max(C_air_new,C_amb_new);

% U_rap = (m_gas/m_air)^0.8;
% 
% NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int)); %da rivedere perché il caso aria-aria non può essere un tubo concentrico

Qid_int = -(T_gas_in - Tamb_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

if T_gas_in<-20

    T_gas_out1 = T_gas_in + Q_int/(cp_air*m_gas);  % kW (K/kW)
    Tamb_out = Tamb_in - Q_int/(cp_amb*m_heater_amb);
    
    Qheat_amb  = Q_int;
    Eps_amb    = epsilon_int;

else 

    T_gas_out1 = T_gas_in;  % kW (K/kW)
    Tamb_out = Tamb_in;
    
    Qheat_amb  = 0;
    Eps_amb    = 0;

end

%% Condizioni di progetto aria-acqua 

T_air_in2 = 0; % °C
% m_air = 3.05; % kg/s

T_water_in = 30;
m_water = 5;

% delta_Ta = 10;

% cp_air = 1.4404; % kJ/ kg K
cp_water = 4.187; % kJ/ kg K

Cmin_ref = min(cp_water*m_water,cp_air*m_air);
Cmax_ref = max(cp_water*m_water,cp_air*m_air);

T_air_out2 = T_water_in - delta_Ta;

Q_rif = m_air*cp_air*(T_air_out2 - T_air_in2);

T_water_out = T_water_in - Q_rif/(cp_water*m_water);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in2 - T_water_out);
DTb = -(T_air_out2 - T_water_in);
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %W/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_rif = UA_rif/Cmin_ref;
% return

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_water_new = cp_water*m_heater_water; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_water_new);
Cmax_new = max(C_air_new,C_water_new);

U_rap = (m_gas/m_air)^0.8;

NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = -(T_gas_out1 - Tw_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

if T_gas_in<Tw_in

    T_gas_out2 = T_gas_out1 + Q_int/(cp_air*m_gas);  % kW (K/kW)
    Tw_out = Tw_in - Q_int/(cp_water*m_heater_water);
    
    Qheat_water  = Q_int;
    Eps_water    = epsilon_int;

else 

    T_gas_out2 = T_gas_out1;  % kW (K/kW)
    Tw_out = Tw_in;
    
    Qheat_water  = 0;
    Eps_water    = 0;

end

%% Condizioni di progetto aria-olio (solare)

T_air_in3 = T_air_out2; % °C
% m_air = 3.2; % kg/s

T_sun_in = 110;
% m_sun = 6;
DTs = 20;
T_sun_out = T_sun_in - DTs;

% delta_Ta = 10;

% cp_air = 1.4404; % kJ/ kg K
cp_sun = 2.195; % kJ/ kg K

T_air_out3 = T_sun_in - delta_Ta;

Q_rif = m_air*cp_air*(T_air_out3 - T_air_in3);

% T_sun_out = T_sun_in - Q_rif/(cp_sun*m_sun);
m_sun = Q_rif/(DTs*cp_sun);

Cmin_ref = min(cp_sun*m_sun,cp_air*m_air);
Cmax_ref = max(cp_sun*m_sun,cp_air*m_air);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in3 - T_sun_out);
DTb = -(T_air_out3 - T_sun_in);
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

Qid_int = -(T_gas_out2 - Ts_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out3 = T_gas_out2 + Q_int/(cp_air*m_gas); % kW (K/kW)
Ts_out = Ts_in - Q_int/(cp_sun*m_heater_sun);

Qheat_sun  = Q_int;
Eps_sun    = epsilon_int;

%% Condizioni di progetto aria-olio 

T_air_in4 = T_air_out3; % °C
% m_air = 3.2; % kg/s

T_oil_in = 230;
% m_oil = 6.4;
DTo = 40;
T_oil_out = T_oil_in - DTo;

% delta_Ta = 10;

% cp_air = 1.4404; % kJ/ kg K
cp_oil = 2.195; % kJ/ kg K

T_air_out4 = 220;

Q_rif = m_air*cp_air*(T_air_out4 - T_air_in4);

% T_oil_out = T_oil_in - Q_rif/(cp_oil*m_oil);
m_oil = Q_rif/(DTs*cp_oil);

Cmin_ref = min(cp_oil*m_oil,cp_air*m_air);
Cmax_ref = max(cp_oil*m_oil,cp_air*m_air);

% METODO DEL DT MEDIO LOGARITMICO
DTa = -(T_air_in4 - T_oil_out);
DTb = -(T_air_out4 - T_oil_in);
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

Qid_int = -(T_gas_out3 - To_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out4 = T_gas_out3 + Q_int/(cp_air*m_gas);  % kW (K/kW)
To_out = To_in - Q_int/(cp_oil*m_heater_oil);

Qheat_oil  = Q_int;
Eps_oil    = epsilon_int;

