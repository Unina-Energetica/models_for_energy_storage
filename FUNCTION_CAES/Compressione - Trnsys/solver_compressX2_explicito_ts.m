function [mc_c1,mc_c2,eta_1,eta_2,Matrix_th,V_Pshaft,massa,Nc_c2] = solver_compressX2_explicito_ts(f_Eta_c1,f_mappa2,f_Eta_c2,beta_1,beta_2,Nc_c1,Mappa_c1_corr,Mappa_c2_corr,parameter,T_oil,T_water,m_oil,m_water)

% date le condizioni di uscita di FSOLVE, calcola il punto di funzionamento
% della coppia di compressori

%% dati ARIA %%

cp_a1 = 1.0426;
cp_a2 = 1.0661;
% cp_a3 = 1.2066;
eps   = (1.4-1)/1.4;

%% delcaration %%

To = parameter(1,1); % bounadary e normalizzation conditions
Ti = parameter(1,2);
po = parameter(1,3);
pi = parameter(1,4);

global graf_e

% main %
[V_mc_C1,V_beta_C1] = single_test_compress(Nc_c1,Mappa_c1_corr);

mc_c1 = interp1(V_beta_C1,V_mc_C1,beta_1,'makima');

if graf_e == 1 
    
    figure(101)
    plot(Mappa_c1_corr(:,:,1),Mappa_c1_corr(:,:,2))
    hold on
    grid on
    plot(V_mc_C1,V_beta_C1,'b--o')
    plot(mc_c1,beta_1,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off 

end

eta_1 = f_Eta_c1(mc_c1,beta_1);

corr = (sqrt(Ti/To))/(pi/po);
m_c1 = mc_c1/corr;

massa = m_c1;

press1 = pi*beta_1;
%% Lavoro di Compressione%%
% il gas si riscalda 

% Tin_c1  = Ti;
% Tout_c1 = Tin_c1*(1 + (beta_1^eps - 1)/eta_1);


%% Intercooler %% #1
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

% cp_cooler   = 4.186; % water
% cp_gas      = cp;  
% T_in_cooler = T_chilling;  % water
% 
% m_gas       = m_c1;
% m_cooler    = m_chill;
% 
% Tin_HE_c1  = Tout_c1;  
% 
% [Tcooler1,Tgas,Qcool1,Eps_HE1] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c1,m_cooler,m_gas);
% 
% Tout_HE_c1 = Tgas; % K

Tin_HE_c1  = Ti*(1 + (beta_1^eps - 1)/eta_1)-273.15; 

[T_gas_out1,To_out1,Qcool_oil1,Eps_oil1,Tout_HE_c1,Tw_out1,Qcool_water1,Eps_water1] = intercooler_centrif_1(T_oil,Tin_HE_c1,T_water,m_oil,m_c1,m_water);

%% Compressore C2 %%

Tin_C2 = Tout_HE_c1+273.15; %ho forzato il raffreddamento per adesso (= Tout_HE_c1)

if Tin_C2 < 0
    Tin_C2 = -Tin_C2;
end

m_c2 =  m_c1;  % rispetto del bilancio di massa tra i due compressori

corr2 = (sqrt(Tin_C2/To))/(press1/po);
mc_c2 = m_c2*corr2;% portata coretta

Nc_c2 = f_mappa2(mc_c2,beta_2);
% Nc_2  = Nc_c2*(sqrt(Tin_C2/To));
% questo è il numero di giri reale ma riferito a una mappa in cui le condizioni d'ingresso sono diverse da quelle della prima mappa,
% quindi non la posso più considerare per le funzioni di interpolazione ottenute dallo scattered

% [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,Mappa_c2_corr);

if graf_e == 1 
    
    figure(102)
    plot(Mappa_c2_corr(:,:,1),Mappa_c2_corr(:,:,2))   
    hold on   
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_c2,beta_2,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
    
    % pause(1)
    
end

eta_2 = f_Eta_c2(mc_c2,beta_2);

P_C1 = (m_c1*cp_a1*Ti*(beta_1^eps-1))/(eta_1);
P_C2 = (m_c2*cp_a2*Tin_C2*(beta_2^eps-1))/(eta_2);

V_Pshaft(1,1) = P_C1;
V_Pshaft(2,1) = P_C2;


%% Intercooler %% #2
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

% cp_cooler   = 4.186; % Water
% cp_gas      = cp;  
% T_in_cooler = T_chilling;
% T_in_gas    = Tin_C2*(1 + (beta_2^eps - 1)/eta_2);
% 
% m_gas       = m_c2;
% m_cooler    = m_chill;
% 
% Tin_HE_c2  = T_in_gas;  
% [Tcooler2,Tgas,Qcool2,Eps_HE2] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c2,m_cooler,m_gas);
% 
% Tout_HE_c2 = Tgas; % K

Tin_HE_c2  = Tin_C2*(1 + (beta_2^eps - 1)/eta_2)-273.15; 

[T_gas_out2,To_out2,Qcool_oil2,Eps_oil2,Tout_HE_c2,Tw_out2,Qcool_water2,Eps_water2] = intercooler_centrif_2(T_oil,Tin_HE_c2,T_water,m_oil,m_c2,m_water);

Matrix_th(1,1) = Tin_HE_c1;  Matrix_th(1,2) = Tout_HE_c1;  Matrix_th(1,3) = T_oil;  Matrix_th(1,4) = To_out1;  Matrix_th(1,5) = T_water;  Matrix_th(1,6) = Tw_out1; Matrix_th(1,7) = Qcool_oil1;  Matrix_th(1,8) = Qcool_water1;  Matrix_th(1,9) = Eps_oil1;  Matrix_th(1,10) = Eps_water1;
Matrix_th(2,1) = Tin_HE_c2;  Matrix_th(2,2) = Tout_HE_c2;  Matrix_th(2,3) = T_oil;  Matrix_th(2,4) = To_out2;  Matrix_th(2,5) = T_water;  Matrix_th(2,6) = Tw_out2; Matrix_th(2,7) = Qcool_oil2;  Matrix_th(2,8) = Qcool_water2;  Matrix_th(2,9) = Eps_oil2;  Matrix_th(2,10) = Eps_water2;

end

