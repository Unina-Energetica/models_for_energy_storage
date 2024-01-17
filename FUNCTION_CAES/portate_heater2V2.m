function Y = portate_heater2V2(mtent,T_gas_in,Ts_in,To_in,m_gas,m_heater_sun)
%function per trovare le portate dei fluidi tramite fsolve

% portate di progetto 
% m_heater_sun   = 5.8;
% m_heater_oil   = 7.7; 

[~,~,~,~,T_gas_out2,~,~,~]=heater2V2(T_gas_in,Ts_in,To_in,m_gas,m_heater_sun,mtent);

Y = 240 - T_gas_out2;

end