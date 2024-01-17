function [Tcooler,Tgas,Qcool,Eps_HE] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

%% declaration %%
T_a = T_in_gas;
T_phi_i = T_in_cooler;
cp_ab = cp_gas; % kj/ kg K
cp_phi_12 = cp_cooler; % kj/ kg K
NTU_int   = 3.5;
%% kernel %%
Cp_ab = cp_ab*m_gas; % kj/ K s -->> kW/K
Cp_phi_12 = cp_phi_12*m_cooler; % kj/ K s -->> kW/K

w_int = min(Cp_ab,Cp_phi_12)/max(Cp_ab,Cp_phi_12); % omega

epsilon_int = (1-exp(-(1-w_int)*NTU_int))/(1-w_int*exp(-(1-w_int)*NTU_int));

Qid_int = (T_a - T_phi_i)*min(Cp_ab,Cp_phi_12); % s -->> kW/K * K  = kW

Q_int = Qid_int*epsilon_int;  % kW/K % calore scambiato nello scambiatore

 T_b = T_a - Q_int/(cp_ab*m_gas);  % kW (K/kW)
 T_phi_u = T_phi_i + Q_int/(cp_phi_12*m_cooler);

 %% output declaration %%

Tcooler  = T_phi_u;
Tgas     = T_b;
Qcool    = Q_int;
Eps_HE   = epsilon_int;



end