function [fn,FVAL,EXITFLAG,OUTPUT] = solver_compressX2_tent(beta_1,beta_2,PShaft,mappa1_corr,mappa2_corr,fun2_C1,fun1_C2,fun2_C2,parameter,tent)
% solutore si entra con P_Available e i fattori di compressore e il
% solutore trova la corretta condizione di funziomaento delle coppia
% compressori
%% declaration
global mappa_C1  mappa_C2  To  Ti  pi  po  power  beta_i1  beta_i2  fun_eta_C1  fun_mappa_C2  fun_eta_C2
  

To = parameter(1,1); % bounadary e normalizzation conditions
Ti = parameter(1,2);
po = parameter(1,3);
pi = parameter(1,4);
%% %% %%                %% %% %%            %% %% %%

mappa_C1     = mappa1_corr;
mappa_C2     = mappa2_corr;

fun_eta_C1   = fun2_C1;
fun_mappa_C2 = fun1_C2; 
fun_eta_C2   = fun2_C2;

power   = PShaft;

beta_i1 = beta_1;
beta_i2 = beta_2;

Y0 = tent;
options = optimset('Display','off');
[fn,FVAL,EXITFLAG,OUTPUT] = fsolve(@Treno_compressori_X2_V3,Y0,options);

end



