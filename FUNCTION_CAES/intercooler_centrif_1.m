function [T_gas_out1,To_out,Qcool_oil,Eps_oil,T_gas_out2,Tw_out,Qcool_water,Eps_water] = intercooler_centrif_1(To_in,T_gas_in,Tw_in,m_cooler_oil,m_gas,m_cooler_water)
%% Intercooler centrifugo 1 %% 

 % = 175;
 % = 650;
 % = 5.8;
 % = 2.9;
 % = 8;
 % = 15;

%% Condizioni di progetto Olio

T_air_in = 690-273.15; % K
m_air = 3; % kg/s

T_oil_in = 180;
m_oil = 6;

delta_To = 20;

cp_air = 1.0515; % kJ/ kg K
cp_oil = 2.195; % kJ/ kg K

Cmin_ref = min(cp_oil*m_oil,cp_air*m_air);
Cmax_ref = max(cp_oil*m_oil,cp_air*m_air);

T_air_out = T_oil_in + delta_To;

Q_rif = m_air*cp_air*(T_air_in - T_air_out);

T_oil_out = T_oil_in + Q_rif/(cp_oil*m_oil);

% METODO DEL DT MEDIO LOGARITMICO
DTa = T_air_in - T_oil_out;
DTb = T_air_out - T_oil_in;
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %kW/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_rif = UA_rif/Cmin_ref;

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_oli_new = cp_oil*m_cooler_oil; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_oli_new);
Cmax_new = max(C_air_new,C_oli_new);

if m_gas < 0

    NTU_int = NTU_rif;

else

    U_rap = (m_gas/m_air)^0.8;
    
    NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

end

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = (T_gas_in - To_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out1 = T_gas_in - Q_int/(cp_air*m_gas);  % kW (K/kW)
To_out = To_in + Q_int/(cp_oil*m_cooler_oil);

Qcool_oil = Q_int;
Eps_oil   = epsilon_int;

%% Condizioni di progetto dissipazione acqua

T_air_in2 = T_air_out; % K
% m_air = 2.85; % kg/s

T_w_in = 15;
m_w = 8;

delta_Tw = 10;

% cp_air = 1.0515; % kJ/ kg K
cp_w = 4.187; % kJ/ kg K

Cmin_ref = min(cp_w*m_w,cp_air*m_air);
Cmax_ref = max(cp_w*m_w,cp_air*m_air);

T_air_out2 = T_w_in + delta_Tw;

Q_rif = m_air*cp_air*(T_air_in2 - T_air_out2);

T_w_out = T_w_in + Q_rif/(cp_w*m_w);

% METODO DEL DT MEDIO LOGARITMICO
DTa = T_air_in2 - T_w_out;
DTb = T_air_out2 - T_w_in;
DTml = (DTa-DTb)/log(DTa/DTb);

UA_rif = Q_rif/(DTml); %W/K

% METODO DELL'EFFICIENZA
% Q_id = (T_air_in - T_oil_in)*Cmin;
% Epsilon = Q/Q_id;

NTU_rif = UA_rif/Cmin_ref;

%% kernel %%
C_air_new = cp_air*m_gas; % kj/ K s -->> kW/K
C_w_new = cp_w*m_cooler_water; % kj/ K s -->> kW/K

Cmin_new = min(C_air_new,C_w_new);
Cmax_new = max(C_air_new,C_w_new);

if m_gas < 0

    NTU_int = NTU_rif;

else

    U_rap = (m_gas/m_air)^0.8;
    
    NTU_int = NTU_rif*U_rap*(Cmin_ref/Cmin_new);

end

w_int = Cmin_new/Cmax_new; % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = (T_gas_out1 - Tw_in)*Cmin_new; % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 %% output declaration %%

T_gas_out2 = T_gas_out1 - Q_int/(cp_air*m_gas);  % kW (K/kW)
Tw_out = Tw_in + Q_int/(cp_w*m_cooler_water);

Qcool_water = Q_int;
Eps_water   = epsilon_int;


end