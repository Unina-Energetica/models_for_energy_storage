function [fn,FVAL,EXITFLAG,OUTPUT] = solver_turbine_X3(P_out,beta_1,beta_2,beta_3,mappa1_corr,mappa2_corr,mappa3_corr,fun_eta_t1,fun_eta_t2,fun_eta_t3,parameter,tent)
% solutore si entra con P_out che si vuole ottenere e i fattori di 
% compressore a disposizione e il solutore trova la corretta condizione di 
% funziomaento della coppia di turbine
%% declaration
global mappa_t1  mappa_t2  mappa_t3  To  Ti  pi  po  power  beta_i1  beta_i2  beta_i3  fun_eta_T1  fun_eta_T2  fun_eta_T3
  

To = parameter(1,1); % bounadary e normalizzation conditions
Ti = parameter(1,2);
po = parameter(1,3);
pi = parameter(1,4);
%% %% %%                %% %% %%            %% %% %%

mappa_t1     = mappa1_corr;
mappa_t2     = mappa2_corr;
mappa_t3     = mappa3_corr;

fun_eta_T1   = fun_eta_t1;
fun_eta_T2   = fun_eta_t2;
fun_eta_T3   = fun_eta_t3;

power   = P_out;

beta_i1 = beta_1;
beta_i2 = beta_2;
beta_i3 = beta_3;

Y0 = [tent];
options = optimset('Display','off');
[fn,FVAL,EXITFLAG,OUTPUT] = fsolve(@Treno_turbine_X3,Y0,options);

end



