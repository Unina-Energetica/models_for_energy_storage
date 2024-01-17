function [MyResults,mFileErrorCode]=Call_compressV2(MyInputs)


global graf_c2 graf_c3 graf_e

N = 2; %numero di compressori cetrifughi
beta_pistoni = 2;
Nsol = 30;

% Parametro per scegliere se stampare i grafici oppure no, non plottando la
% simulazione è più veloce
graf_c2 = 0; %Solo se = 1 stampa.
graf_c3 = 0;
graf_e = 0;

load('mappa_c1_enh.mat'); load('mappa_c2_enh.mat'); %load('mappa_c3_enh.mat')
load('M_ch_su_line_c1.mat'); load('M_ch_su_line_c2.mat'); %load('M_ch_su_line_c3.mat')

load('mappa_c1_corr_enh.mat'); load('mappa_c2_corr_enh.mat'); %load('mappa_c3_corr_enh.mat')
load('M_ch_su_line_c1_corr.mat'); load('M_ch_su_line_c2_corr.mat'); %load('M_ch_su_line_c3_corr.mat')

Mappa_c1 = mappa_c1_enh;
Mappa_c2 = mappa_c2_enh;
% Mappa_c3 = mappa_c3_enh;

Mappa_c1_corr = mappa_c1_corr_enh;
Mappa_c2_corr = mappa_c2_corr_enh;
% Mappa_c3_corr = mappa_c3_corr_enh;

MC_1      = Mappa_c1(:,:,1);
Beta_1    = Mappa_c1(:,:,2);
Nc_1      = Mappa_c1(:,:,3);
ETA_1     = Mappa_c1(:,:,4);

MC_c1      = Mappa_c1_corr(:,:,1);
Beta_c1    = Mappa_c1_corr(:,:,2);
Nc_c1      = Mappa_c1_corr(:,:,3);
ETA_c1     = Mappa_c1_corr(:,:,4);

MC_2      = Mappa_c2(:,:,1);
Beta_2    = Mappa_c2(:,:,2);
Nc_2      = Mappa_c2(:,:,3);
ETA_2     = Mappa_c2(:,:,4);

MC_c2      = Mappa_c2_corr(:,:,1);
Beta_c2    = Mappa_c2_corr(:,:,2);
Nc_c2      = Mappa_c2_corr(:,:,3);
ETA_c2     = Mappa_c2_corr(:,:,4);

% MC_3      = Mappa_c3(:,:,1);
% Beta_3    = Mappa_c3(:,:,2);
% Nc_3      = Mappa_c3(:,:,3);
% ETA_3     = Mappa_c3(:,:,4);
% 
% MC_c3      = Mappa_c3_corr(:,:,1);
% Beta_c3    = Mappa_c3_corr(:,:,2);
% Nc_c3      = Mappa_c3_corr(:,:,3);
% ETA_c3     = Mappa_c3_corr(:,:,4);

choke_line_c1 = M_ch_su_line_c1(:,:,2);
choke_line_c2 = M_ch_su_line_c2(:,:,2);
% choke_line_c3 = M_ch_su_line_c3(:,:,2);

choke_line_c1_corr = M_ch_su_line_c1_corr(:,:,2);
choke_line_c2_corr = M_ch_su_line_c2_corr(:,:,2);
% choke_line_c3_corr = M_ch_su_line_c3_corr(:,:,2);

surge_line_c1 = M_ch_su_line_c1(:,:,1);
surge_line_c2 = M_ch_su_line_c2(:,:,1);
% surge_line_c3 = M_ch_su_line_c3(:,:,1);

surge_line_c1_corr = M_ch_su_line_c1_corr(:,:,1);
surge_line_c2_corr = M_ch_su_line_c2_corr(:,:,1);
% surge_line_c3_corr = M_ch_su_line_c3_corr(:,:,1);

%% data delcaration %%

To = 273.15; %Ti = 15+273.15; %K (stesso valore scelto nella scalatura della mappa)
po = 1.013;  %pi = 1.013; %bar
Ti = MyInputs(3); % K
pi = MyInputs(4); % bar

CINC(1,1) = To; % Compressor inlet and normalizzation conditions
CINC(1,2) = Ti;
CINC(1,3) = po;
CINC(1,4) = pi;

% Pshaft = 2e+03; % Pavalilable = Prenew*eta_inv*eta_reg
% dt = 1; %s
% time_span = 0:dt:3600*3; %s

%Creo delle funzioni interpolanti che, dati i valori di beta e n,
%restituisce i valore di portata corrispondente
f_mc1  = scatteredInterpolant(Beta_1(:),Nc_1(:),MC_1(:));
f_Eta1 = scatteredInterpolant(MC_1(:),Beta_1(:),ETA_1(:));

f_mc_c1  = scatteredInterpolant(Beta_c1(:),Nc_c1(:),MC_c1(:));
f_Nc_c1  = scatteredInterpolant(ETA_c1(:),Beta_c1(:),Nc_c1(:));
f_Eta_c1 = scatteredInterpolant(MC_c1(:),Beta_c1(:),ETA_c1(:));

f_Nc2  = scatteredInterpolant(MC_2(:),Beta_2(:),Nc_2(:));
f_Eta2 = scatteredInterpolant(MC_2(:),Beta_2(:),ETA_2(:));

f_Nc_c2  = scatteredInterpolant(MC_c2(:),Beta_c2(:),Nc_c2(:));
f_Eta_c2 = scatteredInterpolant(MC_c2(:),Beta_c2(:),ETA_c2(:));

% f_Nc3  = scatteredInterpolant(MC_3(:),Beta_3(:),Nc_3(:));
% f_Eta3 = scatteredInterpolant(MC_3(:),Beta_3(:),ETA_3(:));
% 
% f_Nc_c3  = scatteredInterpolant(MC_c3(:),Beta_c3(:),Nc_c3(:));
% f_Eta_c3 = scatteredInterpolant(MC_c3(:),Beta_c3(:),ETA_c3(:));

mFileErrorCode = 110;    % After processing inputs

%Costanti Gas
% cp = 1.004;
% cv = 0.717;
% R = 0.28705; %kJ/(kg*K)

Patm = 1.013; %bar
Beta_plant_lim = 200; %Valore di Beta per il quale scelgo di passare da 2 a 3 compressori
Press_max = 350;

V_power_shaft = zeros(3,1);
Matrix_th = zeros(3,10);
eta_3 = 0;

% Ptk = 90*1.013; %bar
% Vtk = 100; %m3
% T_tk_old = 15+273.15; %valore iniziale di temperatura
% M_tk_old = (Vtk*Ptk*10^2)/(R*(T_tk_old)); %moltiplico per 100 per passare ai kPa, perché R è espresso in kJ
% Tfeed = 25+273.15; 
% press_new = Ptk*100; %kPa
% Pshaft0 = Pshaft;

% ngraf = 10;

% FVAL = zeros(size(time_span)); 
% mc_c1 = zeros(size(time_span)); 
% mc_c2 = zeros(size(time_span));
% % mc_c3 = zeros(size(time_span));
% eta_1 = zeros(size(time_span)); 
% eta_2 = zeros(size(time_span)); 
% eta_3 = zeros(size(time_span)); 
% massa = zeros(size(time_span));
% V_Power = zeros(size(time_span));  
% V_beta1 = zeros(size(time_span));  
% V_beta2 = zeros(size(time_span));
% V_beta3 = zeros(size(time_span));
% eta_average = zeros(size(time_span));
% controllo = zeros(size(time_span));
% Temp_TK = zeros(size(time_span)); 
% V_Tfeed = zeros(size(time_span));
% 
% Matrix_th = zeros(3,6,length(time_span)); 
% V_power_shaft = zeros(3,length(time_span)); 

%% Calcolo Funzionamento %%

% tic

% cont = 1;
% DATA(cont,1) = Ptk/Patm;    %Pressione all'istante t, quando è avvenuta l'espansione 
% DATA(cont,2) = 0;           %Flusso di massa
% DATA(cont,3) = 0;           %Pressione all'ingresso della turbina 2
% DATA(cont,4) = 0;           %Temperatura uscita turbina 1 e in ingresso turbina 2
% DATA(cont,5) = 0;           %Pressione all'ingresso della turbiana 3
% DATA(cont,6) = 0;           %Temperatura uscita turbina 2 e in ingresso turbina 3
% DATA(cont,7) = 0;           %Temperatura uscita turbina 3
% DATA(cont,8) = 0;           %Istante di tempo

% i = 1;

Pel_exc    = MyInputs(1); % kW
Press_tk   = MyInputs(2); % bar
T_oil      = MyInputs(5); % °C
T_water    = MyInputs(6); % °C
m_oil      = MyInputs(7)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s; 
m_water    = MyInputs(8)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s;
%    = MyInputs(9);
% T_chilling  = MyInputs(5); % K
% m_chill     = MyInputs(6); % kg/s

Pshaft = Pel_exc*0.98; %per tenere conto dell'efficenza del motore
Beta_plant = Press_tk/Patm;

if Press_tk < Press_max && Pshaft > 1000 
    
    if Beta_plant < Beta_plant_lim

        % N = 2;
        beta_1     = Beta_plant^(1/N);
        beta_2     = beta_1;
        beta_3     = 0;
    
        %% compressori %%
    
        %% soluzione coppia di compressori %%

        eta_tent = 0.75;
        Ntent = f_Nc_c1(eta_tent,beta_1); % Numero di giri corretto di tentativo, lo scelgo cercando di stare sempre nella zona centrale del grafico
        if Ntent > Mappa_c1_corr(1,end,3)
            % disp('correzione Ntent')
            Ntent = Mappa_c1_corr(1,end,3);
            % pause(1)
        elseif Ntent < Mappa_c1_corr(1,1,3)
            % disp('correzione Ntent')
            Ntent = Mappa_c1_corr(1,1,3);
            % pause(1)
        end 
        mFileErrorCode = 141;
        [NC_c1,~,flag,~] = solver_compressX2_tent_ts(beta_1,beta_2,Pshaft,Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,CINC,Ntent,T_oil,T_water,m_oil,m_water);
        mFileErrorCode = 142;
        if flag <=0 
            mFileErrorCode = 143;
            [pass,NC_c1] = solver_etac_bestX2(Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,Nsol,Pshaft,CINC,T_oil,T_water,m_oil,m_water);
            mFileErrorCode = 144;
            if pass == 1 

                flag = 2;

            else

                %neanche così ho trovato la soluzione, quindi mando
                %l'energia in rete

                mc_c1 = 0;
                mc_c2 = 0;
                eta_1 = 0;
                eta_2 = 0;
                eta_3 = 0;
                massa = 0;
                NC_c2 = 0;
                Matrix_th(:,:)   = 0;
                V_power_shaft(:) = 0;
    
                P_toGrid = Pshaft/0.98;

            end

        end

    %% definizione istante i-esimo
        if flag>0
            mFileErrorCode = 145;
            [mc_c1,mc_c2,eta_1,eta_2,Matrix_th(1:2,:),V_power_shaft(1:2),massa,NC_c2] = solver_compressX2_explicito_ts(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,NC_c1,Mappa_c1_corr,Mappa_c2_corr,CINC,T_oil,T_water,m_oil,m_water);
            mFileErrorCode = 146;
            [f] = compressor_dynLIMIT_NW (NC_c2,beta_2,Mappa_c2_corr,mc_c2);  %% limiti del compressore 2 %%
            flag2 = -2*isnan(f);
            mFileErrorCode = 147;
            if  flag2 == -2

                % disp('soluz fuori limite')

                [pass,NC_c1] = solver_etac_bestX2(Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,Nsol,Pshaft,CINC,T_oil,T_water,m_oil,m_water);

                if pass == 1 

                    [mc_c1,mc_c2,eta_1,eta_2,Matrix_th(1:2,:),V_power_shaft(1:2),massa,NC_c2] = solver_compressX2_explicito_ts(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,NC_c1,Mappa_c1_corr,Mappa_c2_corr,CINC,T_oil,T_water,m_oil,m_water);
  
                    P_toGrid = Pshaft/0.98 - sum(V_power_shaft(:))/0.98;

                    if P_toGrid < 1

                        P_toGrid = 0;

                    end

                else

                    %neanche così ho trovato la soluzione 
                    pass = 0;
    
                    mc_c1 = 0;
                    mc_c2 = 0;
                    eta_1 = 0;
                    eta_2 = 0;
                    eta_3 = 0;
                    massa = 0;
                    Matrix_th(:,:)   = 0;
                    V_power_shaft(:) = 0;
    
                    P_toGrid = Pshaft/0.98; 

                end

            else

                pass = 1; %confermo di aver trovato la soluzione 

                P_toGrid = Pshaft/0.98 - sum(V_power_shaft(:))/0.98;

                if P_toGrid < 1

                    P_toGrid = 0;

                end

            end

        end

    else
        mFileErrorCode = 147.5;
        % N = 3;        
        beta_1     = (Beta_plant/beta_pistoni)^(1/N);
        beta_2     = beta_1;
        beta_3     = beta_pistoni;
    
        %% compressori con pistoni %%
    
        % flag = -10;
    
        % while flag<=0
        %% soluzione coppia di compressori %%
    
        eta_tent = 0.75;
        Ntent = f_Nc_c1(eta_tent,beta_1); % Numero di giri corretto di tentativo, lo scelgo cercando di stare sempre nella zona centrale del grafico
        if Ntent > Mappa_c1_corr(1,end,3)
            % disp('correzione Ntent')
            Ntent = Mappa_c1_corr(1,end,3);
            % pause(1)
        elseif Ntent < Mappa_c1_corr(1,1,3)
            % disp('correzione Ntent')
            Ntent = Mappa_c1_corr(1,1,3);
            % pause(1)
        end 
        mFileErrorCode = 148.2;
        [NC_c1,~,flag,~]    = solver_compressX3_pistoni_tent_ts(beta_1,beta_2,beta_3,Pshaft,Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,CINC,Ntent,T_oil,T_water,m_oil,m_water);
    %   [fn,FVAL,EXITFLAG,OUTPUT] = solver_compressX3_pistoni_tent(beta_1,beta_2,beta_3,PShaft,mappa1_corr,  mappa2_corr,  fun2_C1, fun1_C2,fun2_C2,parameter,tent)

        if flag <= 0
            mFileErrorCode = 149;
            [pass,NC_c1] = solver_etac_bestX3(Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,Nsol,Pshaft,CINC,T_oil,T_water,m_oil,m_water);
            mFileErrorCode = 150;        
            if pass == 1 
    
                flag = 2;
    
            else
    
                %neanche così ho trovato la soluzione, quindi mando
                %l'energia in rete
    
                mc_c1 = 0;
                mc_c2 = 0;
                eta_1 = 0;
                eta_2 = 0;
                eta_3 = 0;
                massa = 0;
                Matrix_th(:,:)   = 0;
                V_power_shaft(:) = 0;
    
                P_toGrid = Pshaft/0.98;
        
            end

        end
        %% definizione istante i-esimo
        if flag>0
            mFileErrorCode = 151;
            [mc_c1,mc_c2,eta_1,eta_2,eta_3,Matrix_th(:,:),V_power_shaft(:),massa,NC_c2] = solver_compressX3_pistoni_explicito_ts(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,NC_c1,Mappa_c1_corr,Mappa_c2_corr,CINC,T_oil,T_water,m_oil,m_water);
            mFileErrorCode = 152;
            [f2] = compressor_dynLIMIT_NW (NC_c2,beta_2,Mappa_c2_corr,mc_c2);  %% limiti del compressore 2 %%
            flag2 = -2*isnan(f2);
            mFileErrorCode = 153;
            if flag2 == -2 
    
                % disp('soluz fuori limite')
                [pass,NC_c1] = solver_etac_bestX3(Mappa_c1_corr,Mappa_c2_corr,f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,Nsol,Pshaft,CINC,T_oil,T_water,m_oil,m_water);
                mFileErrorCode = 154;
                if pass == 1 
                
                    [mc_c1,mc_c2,eta_1,eta_2,eta_3,Matrix_th(:,:),V_power_shaft(:),massa,NC_c2] = solver_compressX3_pistoni_explicito_ts(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,NC_c1,Mappa_c1_corr,Mappa_c2_corr,CINC,T_oil,T_water,m_oil,m_water);
            
                    P_toGrid = Pshaft/0.98 - sum(V_power_shaft(:))/0.98;

                    if P_toGrid < 1

                        P_toGrid = 0;

                    end
                                                
                else
                
                    %neanche così ho trovato la soluzione, quindi mando
                    %l'energia in rete
        
                    mc_c1 = 0;
                    mc_c2 = 0;
                    eta_1 = 0;
                    eta_2 = 0;
                    eta_3 = 0;
                    massa = 0;
                    Matrix_th(:,:)   = 0;
                    V_power_shaft(:) = 0;
        
                    P_toGrid = Pshaft/0.98;
            
                end
                  
            else
                    
                pass = 1; %confermo di aver trovato la soluzione 
                        
                P_toGrid = Pshaft/0.98 - sum(V_power_shaft(:))/0.98;
                        
                if P_toGrid < 1
                
                    P_toGrid = 0;
                
                end
                
            end
        
        end
        
    end

else

    %Serbatoio gia pieno oppure non ho potenza in surplus
    pass = 0;

    beta_1 = 0;
    beta_2 = 0;
    beta_3 = 0;
    mc_c1 = 0;
    mc_c2 = 0;
    eta_1 = 0;
    eta_2 = 0;
    eta_3 = 0;
    massa = 0;
    Matrix_th(:,:)   = 0;
    V_power_shaft(:) = 0;
    
    P_toGrid = Pshaft/0.98;
            
end
           
%% output collection %%
mFileErrorCode = 155;
Pshaft_c1 = V_power_shaft(1);
Pshaft_c2 = V_power_shaft(2);
Pshaft_c3 = V_power_shaft(3);
    
Pw_tot = sum(V_power_shaft); 
    
eta_average = (V_power_shaft(1)*eta_1 + V_power_shaft(2)*eta_2 + V_power_shaft(3)*eta_3)/(Pw_tot);

if isnan(eta_average)

    eta_average = 0;

end

if Beta_plant < Beta_plant_lim 
    T_air_2TK = Matrix_th(2,2);
else        
    T_air_2TK = Matrix_th(3,2);       
end

if pass == 1 
    
    T_olio_sc1    = Matrix_th(1,4);
    T_olio_sc2    = Matrix_th(2,4);
    T_acqua_sc1   = Matrix_th(1,6);
    T_acqua_sc2   = Matrix_th(2,6);
    T_acqua_sc3   = Matrix_th(3,6);
    
else
    
    T_olio_sc1    = T_oil;
    T_olio_sc2    = T_oil;
    T_acqua_sc1   = T_water;
    T_acqua_sc2   = T_water;
    T_acqua_sc3   = T_water;
    
end

Tin_aria_sc1  = Matrix_th(1,1);
Tin_aria_sc2  = Matrix_th(2,1);
Tin_aria_sc3  = Matrix_th(3,1);
Tout_aria_sc1  = Matrix_th(1,2);
Tout_aria_sc2  = Matrix_th(2,2);
Tout_aria_sc3  = Matrix_th(3,2);
Q_olio_sc1    = Matrix_th(1,7);
Q_olio_sc2    = Matrix_th(2,7);
Q_acqua_sc1   = Matrix_th(1,8);
Q_acqua_sc2   = Matrix_th(2,8);
Q_acqua_sc3   = Matrix_th(3,8);
eps_olio_sc1  = Matrix_th(1,9);
eps_olio_sc2  = Matrix_th(2,9);
eps_acqua_sc1 = Matrix_th(1,10);
eps_acqua_sc2 = Matrix_th(2,10);
eps_acqua_sc3 = Matrix_th(3,10);

MyResults(1)=  T_air_2TK; % K
MyResults(2)=  eta_1; 
MyResults(3)=  eta_2;
MyResults(4)=  eta_3;
MyResults(5)=  eta_average;
MyResults(6)=  massa; % kg/s
MyResults(7)=  mc_c1; % kg/s
MyResults(8)=  mc_c2; % kg/s
MyResults(9)=  Pshaft_c1; % kW
MyResults(10)= Pshaft_c2; % kW
MyResults(11)= Pshaft_c3; % kW
MyResults(12)= Tin_aria_sc1; % °C
MyResults(13)= Tin_aria_sc2; % °C
MyResults(14)= Tin_aria_sc3; % °C
MyResults(15)= T_olio_sc1; % °C
MyResults(16)= T_olio_sc2; % °C
MyResults(17)= T_acqua_sc1; % °C
MyResults(18)= T_acqua_sc2; % °C
MyResults(19)= T_acqua_sc3; % °C
MyResults(20)= Q_olio_sc1; % kW
MyResults(21)= Q_olio_sc2; % kW
MyResults(22)= Q_acqua_sc1; % kW
MyResults(23)= Q_acqua_sc2; % kW
MyResults(24)= Q_acqua_sc3; % kW
MyResults(25)= eps_olio_sc1; 
MyResults(26)= eps_olio_sc2;
MyResults(27)= eps_acqua_sc1; 
MyResults(28)= eps_acqua_sc2; 
MyResults(29)= eps_acqua_sc3; 
MyResults(30)= P_toGrid; % kW
MyResults(31)= pass;
MyResults(32)= Pw_tot/0.98;
MyResults(33)= Tout_aria_sc1; % °C
MyResults(34)= Tout_aria_sc2; % °C
MyResults(35)= Tout_aria_sc3; % °C
MyResults(36)= beta_1; % °C
MyResults(37)= beta_2; % °C
MyResults(38)= beta_3; % °C


