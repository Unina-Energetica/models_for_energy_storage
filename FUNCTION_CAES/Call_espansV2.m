function [MyResults,mFileErrorCode]=Call_espansV2(MyInputs)

global graf graf_e graf_b

% Parametro per scegliere se stampare i grafici oppure no, non plottando la
% simulazione è più veloce
graf   = 0; %Solo se = 1 stampa.
graf_e = 0; %Grafico del solver esplicito
graf_b = 0;

N = 3; %Numero di Turbine
Nsol = 30;

%% Carico Mappe delle turbine

load('mappa_t1_enh.mat'); load('mappa_t2_enh.mat'); load('mappa_t3_enh.mat')
load('M_ch_su_line_t1.mat'); load('M_ch_su_line_t2.mat'); load('M_ch_su_line_t3.mat')

load('mappa_t1_corr_enh.mat'); load('mappa_t2_corr_enh.mat'); load('mappa_t3_corr_enh.mat')
load('M_ch_su_line_t1_corr.mat'); load('M_ch_su_line_t2_corr.mat'); load('M_ch_su_line_t3_corr.mat')

Mappa_t1 = mappa_t1_enh;
Mappa_t2 = mappa_t2_enh;
Mappa_t3 = mappa_t3_enh;

Mappa_t1_corr = mappa_t1_corr_enh;
Mappa_t2_corr = mappa_t2_corr_enh;
Mappa_t3_corr = mappa_t3_corr_enh;

M_t1      = Mappa_t1(:,:,1);
Beta_t1   = Mappa_t1(:,:,2);
N_t1      = Mappa_t1(:,:,3);
ETA_t1    = Mappa_t1(:,:,4);

Mc_t1      = Mappa_t1_corr(:,:,1);
Beta_c_t1  = Mappa_t1_corr(:,:,2);
Nc_t1      = Mappa_t1_corr(:,:,3);
ETA_c_t1   = Mappa_t1_corr(:,:,4);

M_t2      = Mappa_t2(:,:,1);
Beta_t2   = Mappa_t2(:,:,2);
N_t2      = Mappa_t2(:,:,3);
ETA_t2    = Mappa_t2(:,:,4);

Mc_t2      = Mappa_t2_corr(:,:,1);
Beta_c_t2  = Mappa_t2_corr(:,:,2);
Nc_t2      = Mappa_t2_corr(:,:,3);
ETA_c_t2   = Mappa_t2_corr(:,:,4);

M_t3      = Mappa_t3(:,:,1);
Beta_t3   = Mappa_t3(:,:,2);
N_t3      = Mappa_t3(:,:,3);
ETA_3     = Mappa_t3(:,:,4);

Mc_t3      = Mappa_t3_corr(:,:,1);
Beta_c_t3  = Mappa_t3_corr(:,:,2);
Nc_t3      = Mappa_t3_corr(:,:,3);
ETA_c_t3   = Mappa_t3_corr(:,:,4);

choke_line_t1 = M_ch_su_line_t1(:,:,2);
choke_line_t2 = M_ch_su_line_t2(:,:,2);
choke_line_t3 = M_ch_su_line_t3(:,:,2);

choke_line_t1_corr = M_ch_su_line_t1_corr(:,:,2);
choke_line_t2_corr = M_ch_su_line_t2_corr(:,:,2);
choke_line_t3_corr = M_ch_su_line_t3_corr(:,:,2);

surge_line_t1 = M_ch_su_line_t1(:,:,1);
surge_line_t2 = M_ch_su_line_t2(:,:,1);
surge_line_t3 = M_ch_su_line_t3(:,:,1);

surge_line_t1_corr = M_ch_su_line_t1_corr(:,:,1);
surge_line_t2_corr = M_ch_su_line_t2_corr(:,:,1);
surge_line_t3_corr = M_ch_su_line_t3_corr(:,:,1);

%Creo delle funzioni interpolanti che, dati i valori di beta e n,
%restituisce i valore di portata corrispondente

% f_m_t1   = scatteredInterpolant(Beta_t1(:),N_t1(:),M_t1(:));
% f_N_t1   = scatteredInterpolant(ETA_t1(:),Beta_t1(:),N_t1(:));
% f_Eta_t1 = scatteredInterpolant(M_t1(:),Beta_t1(:),ETA_t1(:));

f_mc_t1    = scatteredInterpolant(Beta_c_t1(:),Nc_t1(:),Mc_t1(:));
f_Nc_t1    = scatteredInterpolant(ETA_c_t1(:),Beta_c_t1(:),Nc_t1(:));
f_Eta_c_t1 = scatteredInterpolant(Mc_t1(:),Beta_c_t1(:),ETA_c_t1(:));

% f_m_t2   = scatteredInterpolant(Beta_t2(:),N_t2(:),M_t2(:));
% f_N_t2   = scatteredInterpolant(M_t2(:),Beta_t2(:),N_t2(:));
% f_Eta_t2 = scatteredInterpolant(M_t2(:),Beta_t2(:),ETA_t2(:));

f_mc_t2    = scatteredInterpolant(Beta_c_t2(:),Nc_t2(:),Mc_t2(:));
f_Nc_t2    = scatteredInterpolant(Mc_t2(:),Beta_c_t2(:),Nc_t2(:));
f_Eta_c_t2 = scatteredInterpolant(Mc_t2(:),Beta_c_t2(:),ETA_c_t2(:));

% f_m_t3   = scatteredInterpolant(Beta_t3(:),N_t3(:),M_t3(:));
% f_N_t3   = scatteredInterpolant(M_t3(:),Beta_t3(:),N_t3(:));
% f_Eta_t3 = scatteredInterpolant(M_t3(:),Beta_t3(:),ETA_t3(:));

f_mc_t3    = scatteredInterpolant(Beta_c_t3(:),Nc_t1(:),Mc_t1(:));
f_Nc_t3    = scatteredInterpolant(Mc_t3(:),Beta_c_t3(:),Nc_t3(:));
f_Eta_c_t3 = scatteredInterpolant(Mc_t3(:),Beta_c_t3(:),ETA_c_t3(:));

mFileErrorCode = 110;    % After processing inputs

%% Data Declaretion

Pload           = MyInputs(1);
Ptk             = MyInputs(2);
T_tk_old        = MyInputs(3);
M_tk_old        = MyInputs(4);
Tamb_in         = MyInputs(5);
Tw_in           = MyInputs(6);
Ts_in           = MyInputs(7);
To_in           = MyInputs(8);
m_heater_amb    = MyInputs(9);
m_heater_water  = MyInputs(10)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s
m_heater_sun1   = MyInputs(11)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s
m_heater_sun2   = MyInputs(12)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s
m_heater_sun3   = MyInputs(13)/3600; %per passare dai kg/hr delle pompe trnsys a kg/s
% m_heater_oil1 = 10;  % Portata di progetto
% m_heater_oil2 = 7.7; % Portata di progetto
% m_heater_oil3 = 7.4; % Portata di progetto

%Costanti Gas (DA RIVEDERE)
cp = 1.004;
cv = 0.717;
R = 0.28705; %kJ/(kg*K)

% Inizializzazione Tank
% Ptk = 350; %bar
Vtk = 100; %m3
% T_tk_old = 100 + 273.15; %valore iniziale di temperatura del tank
% M_tk_old = (Vtk*Ptk*10^2)/(R*(T_tk_old)); %moltiplico per 100 per passare ai kPa, perché R è espresso in kJ

% Uscita
P_ob = Pload/0.95; %kW, potenza obiettivo 
% P_ob0 = P_ob;
Patm = 1.013; %bar

To = 273.15;  %K
po = 1.013;   %bar
Ti = 200+273.15; %K
pi = Ptk; %bar

TINC(1,1) = To; % Turbine inlet and normalizzation conditions
TINC(1,2) = Ti;
TINC(1,3) = po;
TINC(1,4) = pi;

% Preallocamento variabili
% dt = 1; %s
% time_span = 0:dt:3600*3; %s 

% mc_t1 = zeros(size(time_span)); 
% mc_t2 = zeros(size(time_span));
% mc_t3 = zeros(size(time_span));
% eta_1 = zeros(size(time_span)); 
% eta_2 = zeros(size(time_span)); 
% eta_3 = zeros(size(time_span)); 
% massa = zeros(size(time_span));
% V_Power = zeros(size(time_span));  
% V_beta1 = zeros(size(time_span));  
% V_beta2 = zeros(size(time_span));
% V_beta3 = zeros(size(time_span));
% eta_average = zeros(size(time_span));
% M_tk_new = zeros(size(time_span));
% T_tk_new = zeros(size(time_span)); 

Matrix_th = zeros(3,16); 
Tout_t = zeros(1,3);
V_power_shaft = zeros(3,1); 

% DATA  = zeros(length(time_span),12);     

%% Funzionamento

% T_tk_new(1) = T_tk_old;
% M_tk_new(1) = M_tk_old;

sol_alt = 0;
Press_min = 120;

mFileErrorCode = 141;

if Ptk > Press_min && Pload > 500 
    
    T_gas_in = T_tk_old - 273.15;

    if Ptk >=130 && Ptk < 150
    
        Ti = 190;
    
    elseif Ptk < 130
    
        Ti = 180;
    
    elseif Ptk > 350

        Ti = 210;

    else
        
        Ti = 200;

    end

    iter = 1;

    Beta_plant = Ptk/Patm;

    while iter>0

        TINC(1,2) = Ti+273.15;

        if sol_alt == 1 
    
            beta_1  =  Beta_plant^(1/N)*1.15;
            beta_2  =  Beta_plant^(1/N);%sqrt(Beta_plant/beta_1);
            beta_3  =  Beta_plant/(beta_1*beta_2);
    
        else
        
            beta_1     = Beta_plant^(1/N);
            beta_2     = beta_1;
            beta_3     = beta_2;
    
        end
        %% soluzione coppia di Turbine %%
        
        mFileErrorCode = 142;
    
        [pass,Nc_T1] = solver_etat_bestX3R(Mappa_t1_corr,Mappa_t2_corr,Mappa_t3_corr,f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,Nsol,P_ob,TINC,Ts_in,To_in,m_heater_sun2,m_heater_sun3);    
                
        mFileErrorCode = 143;
        
        if pass == 0 && sol_alt == 0

            % disp('soluz non trovata, cambio beta')
            sol_alt = 1;

            beta_1  =  Beta_plant^(1/N)*1.15;
            beta_2  =  Beta_plant^(1/N);%sqrt(Beta_plant/beta_1);
            beta_3  =  Beta_plant/(beta_1*beta_2);
            
            mFileErrorCode = 144;
            
            [pass,Nc_T1] = solver_etat_bestX3R(Mappa_t1_corr,Mappa_t2_corr,Mappa_t3_corr,f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,Nsol,P_ob,TINC,Ts_in,To_in,m_heater_sun2,m_heater_sun3);

            mFileErrorCode = 145;
            
            if pass == 0 

                % non ho trovato la soluzione 
                % flag = 5;
                mc_t1 = 0;
                mc_t2 = 0;
                mc_t3 = 0;
                eta_1 = 0;
                eta_2 = 0;
                eta_3 = 0;
                massa = 0;
                Matrix_th(:,:)   = 0;
                Tout_t(1,:)      = 0;
                V_power_shaft(:) = 0;
                m_heater_oil1 = 0;
                m_heater_oil2 = 0;
                m_heater_oil3 = 0;
                
                break

            end

        elseif pass == 0 && sol_alt == 1

            % non ho trovato la soluzione 
            % flag = 5;
            mc_t1 = 0;
            mc_t2 = 0;
            mc_t3 = 0;
            eta_1 = 0;
            eta_2 = 0;
            eta_3 = 0;
            massa = 0;
            Matrix_th(:,:)   = 0;
            Tout_t(1,:)      = 0;
            V_power_shaft(:) = 0;
            m_heater_oil1 = 0;
            m_heater_oil2 = 0;
            m_heater_oil3 = 0;
            
            break

        end

        %% Soluzione Esplicita %%
        
        mFileErrorCode = 146;
        
        [mc_t1,mc_t2,mc_t3,eta_1,eta_2,eta_3,V_power_shaft(:,1),massa,Nc_T2,Nc_T3,Tout_t(1,:),Matrix_th(:,:),m_heater_oil2,m_heater_oil3] = solver_TrenoX3_explicitoR(f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,Nc_T1,Mappa_t1_corr,Mappa_t2_corr,Mappa_t3_corr,TINC,Ts_in,To_in,m_heater_sun2,m_heater_sun3);
        
        mFileErrorCode = 147;
        
        fun = @(mtent)portate_preheater(mtent,T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,massa,m_heater_amb,m_heater_water,m_heater_sun1,Ptk);
        mtent = 8;
        options = optimset('Display','off');
        [m_heater_oil1,~,EXITFLAG,~] = fsolve(fun,mtent,options);
        
        mFileErrorCode = 148.1;
    
        if EXITFLAG > 0 && m_heater_oil1<10     
        
            [T_gas_out1,Tamb_out,Qheat_amb,Eps_amb,T_gas_out2,Tw_out,Qheat_water,Eps_water,T_gas_out3,Ts_out,Qheat_sun,Eps_sun,T_gas_out4,To_out,Qheat_oil,Eps_oil]=preheater1ts(T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,massa,m_heater_amb,m_heater_water,m_heater_sun1,m_heater_oil1);

            Matrix_th(1,1) = T_gas_out1;  Matrix_th(1,2) = Tamb_out;   Matrix_th(1,3) = Qheat_amb;   Matrix_th(1,4) = Eps_amb;   Matrix_th(1,5) = T_gas_out2;   Matrix_th(1,6) = Tw_out;   Matrix_th(1,7) = Qheat_water;   Matrix_th(1,8) = Eps_water;
            Matrix_th(1,9) = T_gas_out3;  Matrix_th(1,10) = Ts_out;    Matrix_th(1,11) = Qheat_sun;  Matrix_th(1,12) = Eps_sun;  Matrix_th(1,13) = T_gas_out4;  Matrix_th(1,14) = To_out;  Matrix_th(1,15) = Qheat_oil;    Matrix_th(1,16) = Eps_oil;
             
        else
    
            m_heater_oil1   = 10;
        
            [T_gas_out1,Tamb_out,Qheat_amb,Eps_amb,T_gas_out2,Tw_out,Qheat_water,Eps_water,T_gas_out3,Ts_out,Qheat_sun,Eps_sun,T_gas_out4,To_out,Qheat_oil,Eps_oil]=preheater1ts(T_gas_in,Tamb_in,Tw_in,Ts_in,To_in,massa,m_heater_amb,m_heater_water,m_heater_sun1,m_heater_oil1);
        
            Matrix_th(1,1) = T_gas_out1;  Matrix_th(1,2) = Tamb_out;   Matrix_th(1,3) = Qheat_amb;   Matrix_th(1,4) = Eps_amb;   Matrix_th(1,5) = T_gas_out2;   Matrix_th(1,6) = Tw_out;   Matrix_th(1,7) = Qheat_water;   Matrix_th(1,8) = Eps_water;
            Matrix_th(1,9) = T_gas_out3;  Matrix_th(1,10) = Ts_out;    Matrix_th(1,11) = Qheat_sun;  Matrix_th(1,12) = Eps_sun;  Matrix_th(1,13) = T_gas_out4;  Matrix_th(1,14) = To_out;  Matrix_th(1,15) = Qheat_oil;    Matrix_th(1,16) = Eps_oil;
             
        end
    
    
        if abs(Ti-T_gas_out4)>0.05
    
            Ti = T_gas_out4;
        
        else
    
            iter = 0;
    
        end

    end

    % if flag == 5
    % 
    %         break
    % 
    % end

else

    %Serbatoio "scarico"
    pass = 0;

    beta_1 = 0;
    beta_2 = 0;
    beta_3 = 0;
    mc_t1 = 0;
    mc_t2 = 0;
    mc_t3 = 0;
    eta_1 = 0;
    eta_2 = 0;
    eta_3 = 0;
    massa = 0;
    Matrix_th(:,:)   = 0;
    Tout_t(1,:)      = 0;
    V_power_shaft(:) = 0;
    m_heater_oil1 = 0;
    m_heater_oil2 = 0;
    m_heater_oil3 = 0;
            
end
           
%% output collection %%

mFileErrorCode = 149;

Pshaft_t1 = V_power_shaft(1);
Pshaft_t2 = V_power_shaft(2);
Pshaft_t3 = V_power_shaft(3);

mFileErrorCode = 150;
    
Pw_tot = sum(V_power_shaft(:));
mFileErrorCode = 150.1;
P_fromGrid = Pload - Pw_tot*0.95;

mFileErrorCode = 151;
    
eta_average = (V_power_shaft(1)*eta_1 + V_power_shaft(2)*eta_2 + V_power_shaft(3)*eta_3)/Pw_tot;

if isnan(eta_average)

    eta_average = 0;

end

mFileErrorCode = 152;

Tin_t1 = Matrix_th(1,13); % °C
Tin_t2 = Matrix_th(2,13); % °C
Tin_t3 = Matrix_th(3,13); % °C
Tout_t1 = Tout_t(1,1); % °C
Tout_t2 = Tout_t(1,2); % °C
Tout_t3 = Tout_t(1,3); % °C

mFileErrorCode = 153;

if pass == 1 
    
    T_aria_sc1  = Matrix_th(1,2); % °C
    T_acqua_sc1 = Tw_in-5; %Matrix_th(1,6); % °C
    T_sun_sc1   = Matrix_th(1,10); % °C
    T_olio_sc1  = Matrix_th(1,14); % °C
    T_sun_sc2   = Matrix_th(2,10); % °C
    T_olio_sc2  = Matrix_th(2,14); % °C
    T_sun_sc3   = Matrix_th(3,10); % °C
    T_olio_sc3  = Matrix_th(3,14); % °C
    
else

    T_aria_sc1  = Tamb_in; % °C
    T_acqua_sc1 = Tw_in; % °C
    T_sun_sc1   = Ts_in; % °C
    T_olio_sc1  = To_in; % °C
    T_sun_sc2   = Ts_in; % °C
    T_olio_sc2  = To_in; % °C
    T_sun_sc3   = Ts_in; % °C
    T_olio_sc3  = To_in; % °C
    
end


mFileErrorCode = 154;

Q_amb_sc1     = Matrix_th(1,3); % kW
Q_acqua_sc1   = Matrix_th(1,7); % kW
Q_sun_sc1     = Matrix_th(1,11); % kW
Q_olio_sc1    = Matrix_th(1,15); % kW
Q_sun_sc2     = Matrix_th(2,11); % kW
Q_olio_sc2    = Matrix_th(2,15); % kW
Q_sun_sc3     = Matrix_th(3,11); % kW
Q_olio_sc3    = Matrix_th(3,15); % kW
eps_amb_sc1   = Matrix_th(1,4); 
eps_acqua_sc1 = Matrix_th(1,8);
eps_sun_sc1   = Matrix_th(1,12); 
eps_olio_sc1  = Matrix_th(1,16); 
eps_sun_sc2   = Matrix_th(2,12); 
eps_olio_sc2  = Matrix_th(2,16);
eps_sun_sc3   = Matrix_th(3,12); 
eps_olio_sc3  = Matrix_th(3,16);

mFileErrorCode = 155;

MyResults(1)=  massa; % kg/s
MyResults(2)=  beta_1; % kg/s
MyResults(3)=  beta_2; % kg/s
MyResults(4)=  beta_3; % kg/s
MyResults(5)=  eta_average;
MyResults(6)=  eta_1;
MyResults(7)=  eta_2; 
MyResults(8)=  eta_3; 
MyResults(9)=  Pshaft_t1; % kW
MyResults(10)= Pshaft_t2; % kW
MyResults(11)= Pshaft_t3; % kW
MyResults(12)= Pw_tot*0.95; % kW
MyResults(13)= P_fromGrid; % kW
MyResults(14)= Tout_t1; % °C
MyResults(15)= Tout_t2; % °C
MyResults(16)= Tout_t3; % °C
MyResults(17)= T_aria_sc1; % °C
MyResults(18)= T_acqua_sc1; % °C
MyResults(19)= T_sun_sc1; % °C
MyResults(20)= T_olio_sc1; % °C
MyResults(21)= T_sun_sc2; % °C
MyResults(22)= T_olio_sc2; % °C
MyResults(23)= T_sun_sc3; % °C
MyResults(24)= T_olio_sc3; % °C
MyResults(25)= m_heater_amb; % Ritorno delle portate si conserva
MyResults(26)= m_heater_water*3600; % Ritorno delle portate si conserva, devo restituirla in kg/hr
MyResults(27)= m_heater_sun1*3600; % Ritorno delle portate si conserva, devo restituirla in kg/hr
MyResults(28)= m_heater_sun2*3600; % Ritorno delle portate si conserva, devo restituirla in kg/hr
MyResults(29)= m_heater_sun3*3600; % Ritorno delle portate si conserva, devo restituirla in kg/hr
MyResults(30)= m_heater_oil1*3600; % kg/s
MyResults(31)= m_heater_oil2*3600; % kg/s
MyResults(32)= m_heater_oil3*3600; % kg/s
MyResults(33)= Q_amb_sc1; % kW
MyResults(34)= Q_acqua_sc1; % kW
MyResults(35)= Q_sun_sc1; % kW
MyResults(36)= Q_olio_sc1; % kW
MyResults(37)= Q_sun_sc2; % kW
MyResults(38)= Q_olio_sc2; % kW
MyResults(39)= Q_sun_sc3; % kW
MyResults(40)= Q_olio_sc3; % kW
MyResults(41)= eps_amb_sc1; 
MyResults(42)= eps_acqua_sc1;
MyResults(43)= eps_sun_sc1; 
MyResults(44)= eps_olio_sc1; 
MyResults(45)= eps_sun_sc2; 
MyResults(46)= eps_olio_sc2;
MyResults(47)= eps_sun_sc3; 
MyResults(48)= eps_olio_sc3;
MyResults(49)= pass;
MyResults(50)= Tin_t1; % °C
MyResults(51)= Tin_t2; % °C
MyResults(52)= Tin_t3; % °C










