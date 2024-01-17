function Y = Treno_turbine_X3(Ntent)
% nota la mappa e la potenza disponibile all'albero definisce portata
% corretta, portata istantanea e numero di giri attuali.
global mappa_t1  mappa_t2  mappa_t3  To  Ti  pi  po  power  beta_i1  beta_i2  beta_i3  fun_eta_T1  fun_eta_T2  fun_eta_T3  graf_t 
%% dati ARIA %%

cp  = 1.004;
eps = (1.4-1)/1.4;

%% %% %% ----------- %% %% %%
%% Kernel %%
%% Turbina 1 %%
if Ntent > mappa_t1(1,end,3)
    Ntent = mappa_t1(1,end,3);
elseif Ntent < mappa_t1(1,1,3)
    Ntent = mappa_t1(1,1,3);
end 
[V_mc_t1,V_beta_t1] = single_test_compress(Ntent,mappa_t1);

mc_t1 = interp1(V_beta_t1,V_mc_t1,beta_i1,'makima');

if graf_t == 1 
    
    figure(99)
    plot(mappa_t1(:,:,1),mappa_t1(:,:,2))
    hold on
    grid on
    plot(V_mc_t1,V_beta_t1,'b--o')
    plot(mc_t1,beta_i1,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off    

end

eta_t1 = fun_eta_T1(mc_t1,beta_i1);

corr1 = (sqrt(Ti/To))/(pi/po);

m_t1  = mc_t1/corr1;

press1 = pi/beta_i1;
    %% Lavoro di Compressione%%
    % il gas si riscalda 

%     Tin_c1  = Ti;
%     Tout_c1 = Tin_c1*(beta_i1^eps)/eta_c1;

%% Turbina T2 %%

% Tin_t2 = Ti*(1 - (1 - beta_i1^(-eps))*eta_t1);
Tin_t2 = 273.15 + 200; %considero rilancio termico

m_t2 = m_t1;  % rispetto del bilancio di massa tra i due compressori

corr3 = (sqrt(Tin_t2/To))/(press1/po);
mc_t2 = m_t2*corr3;% portata coretta

%     Nc_c2    = fun_mappa_C2(mc_c2,beta_i2); %credo non serva la velocità 
    
%     [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,mappa_C2);

if graf_t == 1 

    figure(999)
    plot(mappa_t2(:,:,1),mappa_t2(:,:,2))   
    hold on   
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_t2,beta_i2,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
   
end

eta_t2   = fun_eta_T2(mc_t2,beta_i2);

press2 = press1/beta_i2;
%% Turbina T3 %%

Tin_t3 = Tin_t2*(1 - (1 - beta_i2^(-eps))*eta_t2);

m_t3 = m_t2;  % rispetto del bilancio di massa tra i due compressori

corr3 = (sqrt(Tin_t3/To))/(press2/po);
mc_t3 = m_t3*corr3;% portata coretta

%     Nc_c2    = fun_mappa_C2(mc_c2,beta_i2); %credo non serva la velocità 
    
%     [V_mc_C2,V_beta_C2] = single_test_compress(Nc_c2,mappa_C2);

if graf_t == 1 

    figure(9999)
    plot(mappa_t3(:,:,1),mappa_t3(:,:,2))   
    hold on   
    grid on
    % plot(V_mc_C2,V_beta_C2,'b--o')
    plot(mc_t3,beta_i3,'k*')
    xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
    title('Mappa Scalata Corretta')
    hold off
   
end

eta_t3   = fun_eta_T3(mc_t3,beta_i3);

Tout_t3 = Tin_t3*(1 - (1 - beta_i3^(-eps))*eta_t3);

P_t1 = (m_t1*cp*Ti*(1 - beta_i1^(-eps)))*(eta_t1);
P_t2 = (m_t2*cp*Tin_t2*(1 - beta_i2^(-eps)))*(eta_t2);
P_t3 = (m_t3*cp*Tin_t3*(1 - beta_i3^(-eps)))*(eta_t3);
% pause        

Y = P_t1 + P_t2 + P_t3 - power;



% condizione di chiusura %
