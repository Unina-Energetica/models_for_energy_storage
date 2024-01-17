% SOFC model
clc; clear; close all;

%input parameters
R=8.314; 
F=96485;

p_anode=1;                    %bar
p_cathode=1;                  %bar
xH2i=1; xO2i=0.21; xN2i=0.79; %frazioni molari
Acell=87.7;                   %cm2 area attiva della cella 
ncell=8;                      %number of cells
nstack=50;                    %number of stacks
Uf=0.2;                       %fuel utilisation factor
Ua=0.25;                      %air utilisation factor
Tin_H2=750;
Tin_air=750;

P_SOFC   = 25*nstack:10:950*nstack;  %W
V_ideale = 6;               %V 
I_SOFC   = P_SOFC/V_ideale; %A


for i=1:length(I_SOFC)
%     if I_SOFC<=0||I_SOFC>=8300
%        W_H2=0;
%        W_air=0;
%        T=0;
%        eta=0;
%        Vstack=0;
%        Pfc=0;
%     else
  
    Istack(i)=I_SOFC(i)./nstack;
    Icell(i)=Istack(i);
    jcell(i)=Icell(i)./Acell.*1000;
    z_punto(i)=(jcell(i).*Acell./1000)./(2*F); %moli di idrogeno che reagiscono mol/s 
    nH2_in(i)=z_punto(i)./Uf; %moli di idrogeno in ingresso mol/s
    nO2_in(i)=z_punto(i)./2./Ua;  %moli di ossigeno valutate con i coeff stechiometrici mol/s
    nN2_in(i)=nO2_in(i).*xN2i./xO2i; %moli di azoto in ingresso

    nH2_out(i)=nH2_in(i)-z_punto(i); %moli di H2 out
    nO2_out(i)=nO2_in(i)-(0.5.*z_punto(i)); %moli di O2 out
    nH2O_out(i)=z_punto(i); %moli di H2O out
    nN2_out(i)=nN2_in(i); %moli di azoto in uscita
    xH2Oo=nH2O_out(i)./(nH2_out(i)+nH2O_out(i)); 
    xH2o=nH2_out(i)./(nH2_out(i)+nH2O_out(i));
    xO2o=nO2_out(i)./(nO2_out(i)+nN2_out(i));
    
    %entalpia
    T0=273.15;
    A_H2O=3.47;    A_H2=3.3249;   A_O2=3.639;    A_N2=3.28;
    B_H2O=0.00145; B_H2=0.000422; B_O2=0.000506; B_N2=0.000593;
    D_H2O=1210;    D_H2=8300;     D_O2=-22700;   D_N2=4000;
    %IN
    hH2_in=R.*(A_H2.*(Tin_H2-T0)+B_H2./2.*((Tin_H2.^2)-(T0.^2))-D_H2.*(1./Tin_H2-1./T0));
    hO2_in=R.*(A_O2.*(Tin_air-T0)+B_O2./2.*((Tin_air.^2)-(T0.^2))-D_O2.*(1./Tin_air-1./T0));
    hN2_in=R.*(A_N2.*(Tin_air-T0)+B_N2./2.*((Tin_air.^2)-(T0.^2))-D_N2.*(1./Tin_air-1./T0));
    %OUT
    hH2_out=@(T) R.*(A_H2.*(T-T0)+B_H2./2.*((T.^2)-(T0.^2))-D_H2.*(1./T-1./T0));
    hO2_out=@(T) R.*(A_O2.*(T-T0)+B_O2./2.*((T.^2)-(T0.^2))-D_O2.*(1./T-1./T0));
    hH2O_out=@(T) R.*(A_H2O.*(T-T0)+B_H2O./2.*((T.^2)-(T0.^2))-D_H2O.*(1./T-1./T0));
    hN2_out=@(T) R.*(A_N2.*(T-T0)+B_N2./2.*((T.^2)-(T0.^2))-D_N2.*(1./T-1./T0));
    %gibbs
    T0=298.15;
    gf=-228572;
    hf=-241818;
    A=A_H2O-A_H2-0.5*A_O2; B=B_H2O-B_H2-0.5*B_O2; C=0; D=D_H2O-D_H2-0.5*D_O2;
    g=@(T) (gf-hf)./(R.*T0)+hf./(R.*T)+(T.*A+1/2*T.^2.*B+1/3.*T.^3.*C-1./T.*D-T0.*A-1/2.*T0.^2.*B-1/3.*T0.^3.*C+1./T0.*D)./T-(1/2*C.*T.^2+B.*T-1/2*D./T.^2+A.*log(T)-1/2*T0.^2*C-T0.*B+1/2.*D./T0.^2-A.*log(T0));
    gcell=@(T) g(T).*R.*T;
    %rendimento teorico nel caso di trasformazione completamente reversibile
    etamax=@(T) gcell(T)./hf;
    % potenziale della cella senza irreversibilità a circuito aperto
    E0=@(T) -gcell(T)./(2*F);
    % potenziale della cella a circuito aperto
    Ecell=@(T) E0(T)+R.*T.*log((xH2i.*p_anode).*(xO2i.*p_cathode).^0.5./(xH2Oo.*p_anode))./(2*F);
    %perdite per attivazione
    g_a=1344000; 
    g_c=205100;
    ea_a=100000;
    ea_c=120000;
    i0a=@(T) g_a.*exp(-ea_a./(R.*T))*1000; %mA cm-2
    i0c=@(T) g_c.*exp(-ea_c./(R.*T))*1000;  %mA cm-2
    Vacta=@(T) R.*T./F.*asinh(jcell(i)./(2.*i0a(T)));
    Vactc=@(T) R.*T./(2*F).*asinh(jcell(i)./(2.*i0c(T)));
    Vact=@(T) Vacta(T)+Vactc(T);
    %perdite per concentrazione
    il=1900;
    Vconc=@(T) -R*T/(2*F)*log((1-jcell(i)./il));
    %perdite ohmiche
    sigma0=333.3;
    ea_el=85634;
    del=0.00125;
    sigma_el=@(T) sigma0.*exp(-ea_el./(R.*T));
    rel=@(T) del./sigma_el(T);
    rstack=0.057;
    rtot=@(T) (rel(T)+rstack); %ohm*cm2
    Vohm=@(T) rtot(T).*(jcell(i)./1000);
    %potenziale della cella
    Vcell=@(T) Ecell(T)-Vact(T)-Vconc(T)-Vohm(T);
    %bilancio
    f1=(nH2_in(i).*hH2_in+nO2_in(i).*hO2_in+nN2_in(i).*hN2_in);
    f2=@(T) -(nH2_out(i).*hH2_out(T)+nH2O_out(i).*hH2O_out(T)+nO2_out(i).*hO2_out(T)+nN2_out(i).*hN2_out(T));
    f3=@(T) -(Vcell(T).*(jcell(i)./1000).*Acell);
    f4=-(z_punto(i).*hf);
    energy=@(T) f1+f2(T)+f3(T)+f4;
    
    %fzero
    T=fzero(energy,1000);
    
%       if T<=873||T>=1123
%           
%       W_H2=0;
%       W_air=0;
%       T=0;
%       eta=0;
%       Vstack=0;
%       Pfc=0;
%       
%       else
          
      %intensità di corrente della cella
      Icell(i)=jcell(i)./1000.*Acell;
      %potenza della cella
      pcell(i)=Icell(i).*Vcell(T);
      %Cell voltage
      V_cell(i)=Vcell(T);
      E_cell(i)=Ecell(T);
      V_act(i)=Vact(T);
      V_conc(i)=Vconc(T);
      V_ohm(i)=Vohm(T);
      %potenza stack
      Pstack(i)=pcell(i).*ncell;
      %potenza fuel cell
      Pfc(i)=Pstack(i).*nstack;
      %potenziale dello stack
      Vstack(i)=Vcell(T).*ncell;

      %rendimento della cella 
      LHV_H2=120;                                                          %MJ/kg
      PM_H2=2.016;                                                         %g/mol
      mH2_in(i)=nH2_in(i).*PM_H2;                                          %g/s
      eta(i)=pcell(i)./(mH2_in(i).*LHV_H2*1000)*100;                       %adim.
      %flow consumption rate
      nair_in(i)=nO2_in(i)+nN2_in(i);                                      %mol/s
      W_air(i)=(nair_in(i).*ncell.*nstack*3600)*(R/100000)*T/p_cathode;    %m3/h
      W_H2(i)=(nH2_in(i).*ncell*nstack*3600)*(R/100000)*T/p_anode;         %m3/h
%       end

RES.W_H2(i) = W_H2(i);
RES.W_air(i) = W_air(i);
RES.T(i) = T;
RES.eta(i) = eta(i);
RES.V(i) = Vstack(i);
RES.Pfc(i) = Pfc(i)/1000;

end


jcell_start=jcell(1);
jcell_end=jcell(end);
I_SOFC_start=I_SOFC(1);
I_SOFC_end=I_SOFC(end);

figure(1)
plot(jcell,RES.T,'LineWidth',2)
grid on
set(gca,'FontSize',20,'FontAngle','Normal','FontName','Times New Roman');
leg = legend ('$U_f = 0.2$');
set(leg,'Interpreter','latex');
set(leg,'Location','NorthWest','FontSize',20,'FontName','Times New Roman','Orientation','horizontal')
xlabel('i_{cell} [mA/cm^2]')
ylabel('T_{cell} [K]')
xlim([jcell_start jcell_end])

figure(2)
plot(jcell,V_cell,'r')
title('Polarization curve')
xlabel('density current [mA/cm2]')
ylabel('cell voltage [V]')
xlim([jcell_start jcell_end])
grid on
legend('V')
hold on

figure(3)
plot(jcell,RES.W_H2)
hold on
plot(jcell,RES.Pfc)
hold on
xlabel('density current [mA/cm2]')
ylabel('SOFC power [kW]')
yyaxis('right')
ylabel('\eta_{SOFC}')
plot(jcell,RES.eta)
xlim([jcell_start jcell_end])
grid on
hold on
legend('V_{H2} [m3/h]','Pel [kW]','\eta_{SOEC}')
legend('Pel [kW]','\eta_{SOEC}')

figure(4)
plot(jcell,V_cell)
hold on
plot(jcell,E_cell)
hold on
plot(jcell,V_conc)
hold on
plot(jcell,V_ohm)
hold on 
plot(jcell,V_act)
xlabel('density current [mA/cm^{2}]')
ylabel('cell voltage [V]')
xlim([jcell_start jcell_end])
grid on
legend('Vcell','Ecell','Vconc','Vohm','Vact')

figure(5)
plot(jcell,Vstack)
xlabel('density current [mA/cm^{2}]')
ylabel('SOFC voltage [V]')
xlim([jcell_start jcell_end])
grid on
