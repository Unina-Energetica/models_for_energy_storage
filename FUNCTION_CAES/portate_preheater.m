function Y = portate_preheater(mtent,T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,m_gas,m_heater_amb,m_heater_water,m_heater_sun)
%function per trovare le portate dei fluidi tramite fsolve

% global  T_tk_old  Tamb_in  Tw_in  Ts_in  To_in  m_gas
% mtent = [m_heater_water,m_heater_sun,m_heater_oil];
% m_heater_amb   = 30; %portate di progetto
% m_heater_water = 5;
% m_heater_sun   = 8;
% m_heater_oil   = 10;

[~,~,~,~,~,~,~,~,~,~,~,~,T_gas_out4,~,~,~]=preheater1ts(T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,m_gas,m_heater_amb,m_heater_water,m_heater_sun,mtent);
% [T_gas_out1,Tw_out,Qcool_water,Eps_water,T_gas_out2,Ts_out,Qcool_sun,Eps_sun,T_gas_out3,To_out,Qcool_oil,Eps_oil] = preheater1(T_tk_old,Tw_in,Ts_in,To_in,m_gas,m_heater_water,m_heater_sun,m_heater_oil);

    
    Y = 200 - T_gas_out4;

end