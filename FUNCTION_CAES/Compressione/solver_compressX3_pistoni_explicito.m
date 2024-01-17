function [mc_c1,mc_c2,eta_1,eta_2,eta_3,Matrix_th,V_Pshaft,massa,Nc_c2] = solver_compressX3_pistoni_explicito(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,Nc_c1,Mappa_c1_corr,Mappa_c2_corr,parameter)

% date le condizioni di uscita di FSOLVE, calcola il punto di funzionamento
% della coppia di compressori

%% dati ARIA %%

cp  = 1.004;
eps = (1.4-1)/1.4;

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
%     pause

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

cp_cooler   = 4.186; % water
cp_gas      = cp;  
T_in_cooler = 20+273.15;  % water

m_gas       = m_c1;
m_cooler    = 15*m_c1;

Tin_HE_c1  = Ti*(1 + (beta_1^eps - 1)/eta_1); 

[Tcooler1,Tgas,Qcool1,Eps_HE1] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c1,m_cooler,m_gas);

Tout_HE_c1 = Tgas; % K


%% Compressore C2 %%

Tin_C2 = Tout_HE_c1; %ho forzato il raffreddamento per adesso (= Tout_HE_c1)

m_c2 =  m_c1;  % rispetto del bilancio di massa tra i due compressori

corr2 = (sqrt(Tin_C2/To))/(press1/po);
mc_c2 = m_c2*corr2;% portata coretta

Nc_c2 = f_Nc_c2(mc_c2,beta_2); %numero di giri corretto
% Nc_2  = Nc_c2*(sqrt(Tin_C2/To)); %numero di giri reale, stesso discorso
% in solver esplicito X2

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

end

eta_2 = f_Eta_c2(mc_c2,beta_2);

press2 = press1*beta_2;
%% Intercooler %% #2
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

cp_cooler   = 4.186; % Water
cp_gas      = cp;  
T_in_cooler = 20+273.15;
Tin_HE_c2    = Tin_C2*(1 + (beta_2^eps - 1)/eta_2);

m_gas       = m_c2;
m_cooler    = 15*m_c2;
 
[Tcooler2,Tgas,Qcool2,Eps_HE2] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c2,m_cooler,m_gas);

Tout_HE_c2 = Tgas; % K

%% Compressore C3 pistoni%%

Tin_C3 = Tout_HE_c2; %ho forzato il raffreddamento per adesso (= Tout_HE_c1)

m_c3 =  m_c2;  % rispetto del bilancio di massa tra i due compressori

% corr3 = (sqrt(Tin_C3/To))/(press2/po);
% mc_c3 = m_c3*corr3;% portata coretta

% Nc_c3 = f_Nc_c3(mc_c3,beta_3); %numero di giri corretto
% Nc_3  = Nc_c3*(sqrt(Tin_C2/To)); %numero di giri reale

% if graf_e == 1 
% 
%     figure(103)
%     plot(Mappa_c3_corr(:,:,1),Mappa_c3_corr(:,:,2))
%     hold on
%     grid on
%     % plot(V_mc_C3,V_beta_C3,'b--o')
%     plot(mc_c3,beta_3,'k*')
%     xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
%     title('Mappa Scalata Corretta')
%     hold off 
% 
% %     pause
% 
% end

% eta_3 = f_Eta_c3(mc_c3,beta_3);

eta_3 = 0.82;

P_C1 = (m_c1*cp*Ti*(beta_1^eps-1))/(eta_1);
P_C2 = (m_c2*cp*Tin_C2*(beta_2^eps-1))/(eta_2);
P_C3 = (m_c3*cp*Tin_C3*(beta_3^eps-1))/(eta_3);

V_Pshaft(1,1) = P_C1;
V_Pshaft(2,1) = P_C2;
V_Pshaft(3,1) = P_C3;



%% Intercooler %% #3
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 
% intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,T_in_gas,m_cooler,m_gas)

cp_cooler   = 4.186; % Water
cp_gas      = cp;  
T_in_cooler = 20+273.15;
Tin_HE_c3    = Tin_C3*(1 + (beta_3^eps - 1)/eta_3);

m_gas       = m_c3;
m_cooler    = 15*m_c3;
 
[Tcooler3,Tgas,Qcool3,Eps_HE3] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c3,m_cooler,m_gas);

Tout_HE_c3 = Tgas; % K

Matrix_th(1,1) = Tin_HE_c1;  Matrix_th(1,2) = Tout_HE_c1; Matrix_th(1,3) = T_in_cooler;  Matrix_th(1,4) = Tcooler1; Matrix_th(1,5) = Qcool1; Matrix_th(1,6) = Eps_HE1;
Matrix_th(2,1) = Tin_HE_c2;  Matrix_th(2,2) = Tout_HE_c2; Matrix_th(2,3) = T_in_cooler;  Matrix_th(2,4) = Tcooler2; Matrix_th(2,5) = Qcool2; Matrix_th(2,6) = Eps_HE2;
Matrix_th(3,1) = Tin_HE_c3;  Matrix_th(3,2) = Tout_HE_c3; Matrix_th(3,3) = T_in_cooler;  Matrix_th(3,4) = Tcooler3; Matrix_th(3,5) = Qcool3; Matrix_th(2,6) = Eps_HE3;


%% TK %%