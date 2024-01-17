function Y = Treno_compressori_X3_pistoni(Ntent)
% nota la mappa e la potenza disponibile all'albero definisce portata
% corretta, portata istantanea e numero di giri attuali.
global mappa_C1  mappa_C2  To  Ti  pi  po  power  beta_i1  beta_i2  beta_i3  fun_eta_C1  fun_mappa_C2  fun_eta_C2  graf_c3
         
%% dati ARIA %%

cp  = 1.004;
eps = (1.4-1)/1.4;

%% %% %% ----------- %% %% %%
%% Compressore 1 %%

[V_mc_C1,V_beta_C1] = single_test_compress(Ntent,mappa_C1);

mc_c1 = interp1(V_beta_C1,V_mc_C1,beta_i1,'makima');

if graf_c3 == 1 

    figure(99)
    plot(mappa_C1(:,:,1),mappa_C1(:,:,2))
    hold on
    grid on
    plot(V_mc_C1,V_beta_C1,'b--o')
    plot(mc_c1,beta_i1,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off

end

eta_c1 = fun_eta_C1(mc_c1,beta_i1);

corr = (sqrt(Ti/To))/(pi/po);

m_c1  = mc_c1/corr;

press1 = pi*beta_i1;
%% Lavoro di Compressione%%
% il gas si riscalda 

% Tin_c1  = Ti;
% Tout_c1 = Tin_c1*(1 + (beta_1^eps - 1)/eta_1);

%% Intercooler %%
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 


cp_cooler   = 4.186; %
cp_gas      = cp;  
T_in_cooler = 20+273.15;
Tin_HE_c1  = Ti*(1 + (beta_i1^eps - 1)/eta_c1); 

m_gas       = m_c1;
m_cooler    = 15*m_c1;

[~,Tgas,~,~] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c1,m_cooler,m_gas);

Tout_HE_c1 = Tgas; % K

%% Compressore C2 %%

Tin_C2 = Tout_HE_c1;

if Tin_C2 < 0
    Tin_C2 = -Tin_C2;
end

m_c2 = m_c1;  % rispetto del bilancio di massa tra i due compressori

corr2 = (sqrt(Tin_C2/To))/(press1/po);
mc_c2 = m_c2*corr2;% portata coretta

% Nc_c2    = fun_mappa_C2(mc_c2,beta_i2) %credo non serva la velocità     
% [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,mappa_C2);

if graf_c3 == 1 

    figure(999)
    plot(mappa_C2(:,:,1),mappa_C2(:,:,2))
    hold on
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_c2,beta_i2,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
        
end

eta_c2   = fun_eta_C2(mc_c2,beta_i2);
        
press2 = press1*beta_i2;
%% Lavoro di Compressione%%
% il gas si riscalda 

% Tout_c2 = Tin_C2*(1 + (beta_2^eps - 1)/eta_2);

%% Intercooler %%
% raffreddamento intercooling 
% (attualmente scambio termico semplificato  mettere metodo NTU-eps 

cp_cooler   = 4.186; %
cp_gas      = cp;  
T_in_cooler = 20+273.15;
Tin_HE_c2  = Tin_C2*(1 + (beta_i2^eps - 1)/eta_c2); 

m_gas       = m_c2;
m_cooler    = 15*m_c2;

[~,Tgas,~,~] = intercooler_compressor(cp_cooler,cp_gas,T_in_cooler,Tin_HE_c2,m_cooler,m_gas);

%Per adesso considero intercooler semplificato
 
Tout_HE_c2 = Tgas; % K

%% Compressore C3 pistoni %%

Tin_C3 = Tout_HE_c2;

m_c3 = m_c2;  % rispetto del bilancio di massa tra i due compressori

% corr3 = (sqrt(Tin_C3/To))/(press2/po);
% mc_c3 = m_c3*corr3;% portata coretta

% Nc_c3    = fun_mappa_C3(mc_c3,beta_i3); %credo non serva la velocità 
        
% [V_mc_C3,V_beta_C3] = single_test_compress(Nc_c3,mappa_C3);

% if graf == 1 
% 
%     figure(9999)
%     plot(mappa_C3(:,:,1),mappa_C3(:,:,2))
%     hold on
%     grid on
%     % plot(V_mc_C3,V_beta_C3,'b--o')
%     plot(mc_c3,beta_i3,'k*')
%     xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
%     title('Mappa Scalata')
%     hold off
% 
% end

% eta_c3   = fun_eta_C3(mc_c3,beta_i3);

eta_c3 = 0.82; %da rivedere
            
P_C1 = (m_c1*cp*Ti*(beta_i1^eps-1))/(eta_c1);
P_C2 = (m_c2*cp*Tin_C2*(beta_i2^eps-1))/(eta_c2);
P_C3 = (m_c3*cp*Tin_C3*(beta_i3^eps-1))/(eta_c3);
            
% pause
            
Y = P_C1 + P_C2 + P_C3 - power;


% condizione di chiusura %


