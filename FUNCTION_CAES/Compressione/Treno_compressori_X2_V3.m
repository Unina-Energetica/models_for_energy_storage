function Y = Treno_compressori_X2_V3(Ntent)
% nota la mappa e la potenza disponibile all'albero definisce portata
% corretta, portata istantanea e numero di giri attuali.
global mappa_C1  mappa_C2  To  Ti  pi  po  power  beta_i1  beta_i2  fun_eta_C1  fun_mappa_C2  fun_eta_C2  graf_c2
%% dati ARIA %%

cp  = 1.004;
eps = (1.4-1)/1.4;

%% %% %% ----------- %% %% %%
%% Kernel %%
%% Compressore 1 %%

if Ntent > mappa_C1(1,end,3)
    % disp('correzione Ntent')
    Ntent = mappa_C1(1,end,3);
elseif Ntent < mappa_C1(1,1,3)
    % disp('correzione Ntent');
    Ntent = mappa_C1(1,1,3);
end 

[V_mc_C1,V_beta_C1] = single_test_compress(Ntent,mappa_C1); % Se sceglie un numero di giri troppo basso estrapola un portata inesistente 

mc_c1 = interp1(V_beta_C1,V_mc_C1,beta_i1,'makima');

if graf_c2 == 1 
    
    figure(88)
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

%     Tin_c1  = Ti;
%     Tout_c1 = Tin_c1*(1 + (beta_1^eps - 1)/eta_1);

%% Intercooler %% #1
% raffreddamento intercooling 
% attualmente scambio termico semplificato metodo NTU-eps 

cp_cooler   = 4.186; % water
cp_gas      = cp;  
T_in_cooler = 20+273.15;  % water

m_gas       = m_c1;
m_cooler    = 15*m_c1;

Tin_HE_c1  = Ti*(1 + (beta_i1^eps - 1)/eta_c1);  

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

if graf_c2 == 1 

    figure(888)
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

P_C1 = (m_c1*cp*Ti*(beta_i1^eps-1))/(eta_c1);
P_C2 = (m_c2*cp*Tin_C2*(beta_i2^eps-1))/(eta_c2);

% P_C1
% P_C2
% Ptot = P_C1 + P_C2
% pause(0.5)
        
Y = P_C1 + P_C2 - power;



% condizione di chiusura %
