function [mc_t1,mc_t2,mc_t3,eta_1,eta_2,eta_3,V_Pshaft,massa,Nc_t2,Nc_t3,Tout_t] = solver_TrenoX3_explicito(f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,Nc_t1,Mappa_t1_corr,Mappa_t2_corr,Mappa_t3_corr,parameter) 
% ,Matrix_th
% date le condizioni di uscita di FSOLVE, calcola il punto di funzionamento
% della coppia di compressori

%% dati ARIA %%

cp  = 1.004;
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
%% Lavoro di Espansione%% 

Tin_11  = Ti;
Tout_t(1,1) = Ti*(1 - (1 - beta_1^(-eps))*eta_1);

%% Intercooler %% #1
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

% cp_cooler   = 4.186; % water
% cp_gas      = cp;  
% T_in_cooler = 20+273.15;  % water

% m_gas       = m_t1;
% m_cooler    = 15*m_t1;

% Tin_HE_c1  = Tout_c1;  

% [Tcooler1,Tgas,Qcool1,Eps_HE1] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c1,m_cooler,m_gas);

% Tout_HE_c1 = Tgas; % K


%% Turbina T2 %%

% Tin_t2 = Ti*(1 - (1 - beta_1^(-eps))*eta_1); 
Tin_t2 = 273.15 + 200; %considero rilancio termico

m_t2 =  m_t1;  % rispetto del bilancio di massa tra i due compressori

corr2 = (sqrt(Tin_t2/To))/(press1/po);
mc_t2 = m_t2*corr2;% portata coretta

Nc_t2 = f_Nc_t2(mc_t2,beta_2);
% Nc_2  = Nc_c2*(sqrt(Tin_C2/To));
% questo è il numero di giri reale ma riferito a una mappa in cui le condizioni d'ingresso sono diverse da quelle della prima mappa,
% quindi non la posso più considerare per le funzioni di interpolazione ottenute dallo scattered

% [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,Mappa_c2_corr);

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
%% Turbina T3 %%

Tin_t3 = Tin_t2*(1 - (1 - beta_2^(-eps))*eta_2);
% Tin_t3 = 273.15 + 200; %considero rilancio termico
Tout_t(1,2) = Tin_t3;

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

Tout_t(1,3) = Tin_t3*(1 - (1 - beta_3^(-eps))*eta_3);

P_t1 = (m_t1*cp*Ti*(1 - beta_1^(-eps)))*(eta_1);
P_t2 = (m_t2*cp*Tin_t2*(1 - beta_2^(-eps)))*(eta_2);
P_t3 = (m_t3*cp*Tin_t3*(1 - beta_3^(-eps)))*(eta_3);

V_Pshaft(1,1) = P_t1;
V_Pshaft(2,1) = P_t2;
V_Pshaft(3,1) = P_t3;

%% Intercooler %% #2
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

% cp_cooler   = 4.186; % Water
% cp_gas      = cp;  
% T_in_cooler = 20+273.15;
% T_in_gas    = Tin_T2*(beta_2^eps)/eta_2;
% m_gas       = m_c2;
% m_cooler    = 15*m_c2;

% Tin_HE_c2  = T_in_gas;  
% [Tcooler2,Tgas,Qcool2,Eps_HE2] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c2,m_cooler,m_gas);


% Tout_HE_c2 = Tgas; % K

% Matrix_th(1,1) = Tin_HE_c1;  Matrix_th(1,2) = Tout_HE_c1; Matrix_th(1,3) = T_in_cooler;  Matrix_th(1,4) = Tcooler1; Matrix_th(1,5) = Qcool1; Matrix_th(1,6) = Eps_HE1;
% Matrix_th(2,1) = Tin_HE_c2;  Matrix_th(2,2) = Tout_HE_c2; Matrix_th(2,3) = T_in_cooler;  Matrix_th(2,4) = Tcooler2; Matrix_th(2,5) = Qcool2; Matrix_th(2,6) = Eps_HE2;



%% TK %%