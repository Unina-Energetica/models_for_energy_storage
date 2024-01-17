function [mc_t1,mc_t2,mc_t3,eta_1,eta_2,eta_3,V_Pshaft,massa,Nc_t2,Nc_t3,Tout_t,Matrix_th,m_heater_oil2,m_heater_oil3] = solver_TrenoX3_explicitoRV2(f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,Nc_t1,Mappa_t1_corr,Mappa_t2_corr,Mappa_t3_corr,parameter,Ts_in,To_in,m_heater_sun2,m_heater_sun3) 
%% dati ARIA %%

cp_t1  = 1.0953;
cp_t2  = 1.0298;
cp_t3  = 1.0169;
eps = (1.4-1)/1.4;

%% delcaration %%

To = parameter(1,1); % Inlet and normalizzation conditions
Ti = parameter(1,2);
po = parameter(1,3);
pi = parameter(1,4);

global graf_e

% main %
if Nc_t1 > Mappa_t1_corr(1,end,3)
    Nc_t1 = Mappa_t1_corr(1,end,3);
elseif Nc_t1 < Mappa_t1_corr(1,1,3)
    Nc_t1 = Mappa_t1_corr(1,1,3);
end 
[V_mc_t1,V_beta_t1] = single_test_compress(Nc_t1,Mappa_t1_corr);

mc_t1 = interp1(V_beta_t1,V_mc_t1,beta_1,'makima');

if graf_e == 1 
    
    figure(201)
    plot(Mappa_t1_corr(:,:,1),Mappa_t1_corr(:,:,2))
    hold on
    grid on
    plot(V_mc_t1,V_beta_t1,'b--o')
    plot(mc_t1,beta_1,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off 

end

eta_1 = f_Eta_c_t1(mc_t1,beta_1);

corr = (sqrt(Ti/To))/(pi/po);

m_t1 = mc_t1/corr;
massa = m_t1;

press1 = pi/beta_1;

Tout_t(1,1) = Ti*(1 - (1 - beta_1^(-eps))*eta_1) - 273.15;

%% Heater #2 %%

fun = @(mtent)portate_heater2V2(mtent,Tout_t(1,1),Ts_in,To_in,massa,m_heater_sun2);
mtent = 5;
options = optimset('Display','off');
[m_heater_oil2,~,EXITFLAG,~] = fsolve(fun,mtent,options);

if EXITFLAG > 0 && m_heater_oil2<6     
            
    [T_gas_out1,Ts_out,Qheat_sun,Eps_sun,T_gas_out2,To_out,Qheat_oil,Eps_oil]=heater2V2(Tout_t(1,1),Ts_in,To_in,massa,m_heater_sun2,m_heater_oil2);

  % Matrix_th(2,1) = 0;  Matrix_th(2,2) = 0;   Matrix_th(2,3) = 0;   Matrix_th(2,4) = 0;   Matrix_th(2,5) = 0;   Matrix_th(2,6) = 0;   Matrix_th(2,7) = 0;   Matrix_th(2,8) = 0;
    Matrix_th(2,9) = T_gas_out1;  Matrix_th(2,10) = Ts_out;    Matrix_th(2,11) = Qheat_sun;  Matrix_th(2,12) = Eps_sun;  Matrix_th(2,13) = T_gas_out2;  Matrix_th(2,14) = To_out;  Matrix_th(2,15) = Qheat_oil;    Matrix_th(2,16) = Eps_oil;
     
else

    m_heater_oil2   = 6; % Portata di progetto

    [T_gas_out1,Ts_out,Qheat_sun,Eps_sun,T_gas_out2,To_out,Qheat_oil,Eps_oil]=heater2V2(Tout_t(1,1),Ts_in,To_in,massa,m_heater_sun2,m_heater_oil2);

  % Matrix_th(2,1) = 0;  Matrix_th(2,2) = 0;   Matrix_th(2,3) = 0;   Matrix_th(2,4) = 0;   Matrix_th(2,5) = 0;   Matrix_th(2,6) = 0;   Matrix_th(2,7) = 0;   Matrix_th(2,8) = 0;
    Matrix_th(2,9) = T_gas_out1;  Matrix_th(2,10) = Ts_out;    Matrix_th(2,11) = Qheat_sun;  Matrix_th(2,12) = Eps_sun;  Matrix_th(2,13) = T_gas_out2;  Matrix_th(2,14) = To_out;  Matrix_th(2,15) = Qheat_oil;    Matrix_th(2,16) = Eps_oil;

end

%% Turbina T2 %%

Tin_t2 = T_gas_out2+273.15; %considero rilancio termico

m_t2 =  m_t1;  % rispetto del bilancio di massa tra i due compressori

corr2 = (sqrt(Tin_t2/To))/(press1/po);
mc_t2 = m_t2*corr2;% portata coretta

Nc_t2 = f_Nc_t2(mc_t2,beta_2);

if graf_e == 1 
    
    figure(202)
    plot(Mappa_t2_corr(:,:,1),Mappa_t2_corr(:,:,2))   
    hold on   
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_t2,beta_2,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
    
end

eta_2 = f_Eta_c_t2(mc_t2,beta_2);

press2 = press1/beta_2;

Tout_t(1,2) = Tin_t2*(1 - (1 - beta_2^(-eps))*eta_2) - 273.15;

%% Heater #3 %%

fun = @(mtent)portate_heater3V2(mtent,Tout_t(1,2),Ts_in,To_in,massa,m_heater_sun3);
mtent = 5;
options = optimset('Display','off');
[m_heater_oil3,~,EXITFLAG,~] = fsolve(fun,mtent,options);

if EXITFLAG > 0 && m_heater_oil3<10    
            
    [T_gas_out1,Ts_out,Qheat_sun,Eps_sun,T_gas_out2,To_out,Qheat_oil,Eps_oil]=heater3V2(Tout_t(1,2),Ts_in,To_in,massa,m_heater_sun3,m_heater_oil3);

  % Matrix_th(3,1) = 0;  Matrix_th(3,2) = 0;   Matrix_th(3,3) = 0;   Matrix_th(3,4) = 0;   Matrix_th(3,5) = 0;   Matrix_th(3,6) = 0;   Matrix_th(3,7) = 0;   Matrix_th(3,8) = 0;
    Matrix_th(3,9) = T_gas_out1;  Matrix_th(3,10) = Ts_out;    Matrix_th(3,11) = Qheat_sun;  Matrix_th(3,12) = Eps_sun;  Matrix_th(3,13) = T_gas_out2;  Matrix_th(3,14) = To_out;  Matrix_th(3,15) = Qheat_oil;    Matrix_th(3,16) = Eps_oil;
     
else

    m_heater_oil3   = 10; % Portata di progetto

    [T_gas_out1,Ts_out,Qheat_sun,Eps_sun,T_gas_out2,To_out,Qheat_oil,Eps_oil]=heater3V2(Tout_t(1,2),Ts_in,To_in,massa,m_heater_sun3,m_heater_oil3);

  % Matrix_th(3,1) = 0;  Matrix_th(3,2) = 0;   Matrix_th(3,3) = 0;   Matrix_th(3,4) = 0;   Matrix_th(3,5) = 0;   Matrix_th(3,6) = 0;   Matrix_th(3,7) = 0;   Matrix_th(3,8) = 0;
    Matrix_th(3,9) = T_gas_out1;  Matrix_th(3,10) = Ts_out;    Matrix_th(3,11) = Qheat_sun;  Matrix_th(3,12) = Eps_sun;  Matrix_th(3,13) = T_gas_out2;  Matrix_th(3,14) = To_out;  Matrix_th(3,15) = Qheat_oil;    Matrix_th(3,16) = Eps_oil;
 
end

%% Turbina T3 %%

% Tin_t3 = Tin_t2*(1 - (1 - beta_2^(-eps))*eta_2);
Tin_t3 = T_gas_out2+273.15; %considero rilancio termico

m_t3 = m_t2;  % rispetto del bilancio di massa tra i due compressori

corr3 = (sqrt(Tin_t3/To))/(press2/po);
mc_t3 = m_t3*corr3;% portata coretta

Nc_t3 = f_Nc_t3(mc_t3,beta_3);
%     Nc_c2    = fun_mappa_C2(mc_c2,beta_i2); %credo non serva la velocità 
    
%     [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,mappa_C2);

if graf_e == 1 

    figure(203)
    plot(Mappa_t3_corr(:,:,1),Mappa_t3_corr(:,:,2))   
    hold on   
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_t3,beta_3,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
    
    pause(1)
   
end

eta_3   = f_Eta_c_t3(mc_t3,beta_3);

Tout_t(1,3) = Tin_t3*(1 - (1 - beta_3^(-eps))*eta_3) - 273.15;

P_t1 = (m_t1*cp_t1*Ti*(1 - beta_1^(-eps)))*(eta_1);
P_t2 = (m_t2*cp_t2*Tin_t2*(1 - beta_2^(-eps)))*(eta_2);
P_t3 = (m_t3*cp_t3*Tin_t3*(1 - beta_3^(-eps)))*(eta_3);

V_Pshaft(1,1) = P_t1;
V_Pshaft(2,1) = P_t2;
V_Pshaft(3,1) = P_t3;

end

