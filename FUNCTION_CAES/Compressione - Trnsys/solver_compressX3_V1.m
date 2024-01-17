function [fn,FVAL,EXITFLAG,OUTPUT] = solver_trenoX3_compressV1(beta_1,beta_2,beta_3,PShaft,mappa1_corr,mappa2_corr,mappa3_corr,fun2_C1,fun1_C2,fun2_C2,fun1_C3,fun2_C3,parameter,tent)
% solutore si entra con P_Available e i fattori di compressore e il
% solutore trova la corretta condizione di funziomaento delle coppia
% compressori
%% declaration
global mappa_C1  mappa_C2  mappa_C3  To  Ti  pi  po  power  beta_i1  beta_i2  beta_i3  fun_eta_C1  fun_mappa_C2  fun_eta_C2  fun_mappa_C3  fun_eta_C3  graf
  

To = parameter(1,1); % bounadary e normalizzation conditions
Ti = parameter(1,2);
po = parameter(1,3);
pi = parameter(1,4);
%% %% %%                %% %% %%            %% %% %%

mappa_C1     = mappa1_corr;
mappa_C2     = mappa2_corr;
mappa_C3     = mappa3_corr;

fun_eta_C1   = fun2_C1;
fun_mappa_C2 = fun1_C2; 
fun_eta_C2   = fun2_C2;
fun_mappa_C3 = fun1_C3; 
fun_eta_C3   = fun2_C3;

power   = PShaft;

beta_i1 = beta_1;
beta_i2 = beta_2;
beta_i3 = beta_3;

Y0 = [tent];
[fn,FVAL,EXITFLAG,OUTPUT] = fsolve(@Treno_compressori_X3_V1,Y0);

end



