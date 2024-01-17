function [pass,NC_c1] = solver_etac_bestX3(mappa_c1,mappa_c2,f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,Nsol,P_disp,CINC,T_oil,T_water,m_oil,m_water)

global graf_bc3 

j=0;
for y = 1:length(mappa_c1(1,:,3))

    if beta_1>min(mappa_c1(:,y,2)) && beta_1<max(mappa_c1(:,y,2))
    
        j = j + 1;
        if j == 1
            minN_aly = mappa_c1(1,y,3);
            maxN_aly = mappa_c1(1,y,3);
        end
        
        if mappa_c1(1,y,3)>maxN_aly
            maxN_aly = mappa_c1(1,y,3);
        elseif mappa_c1(1,y,3)<minN_aly
            minN_aly = mappa_c1(1,y,3);
        end
        % N_aly(j) = mappa_c1(1,y,3);  
        % beta_test(j) = beta_1;
    
    end

end

% minN_aly = min(N_aly);
% maxN_aly = max(N_aly);
% clear N_aly

N_aly = linspace(minN_aly,maxN_aly,Nsol);

for y = 1:length(N_aly)
    
    [mc_c1test(y),mc_c2test(y),eta_b1(y),eta_b2(y),eta_b3(y),~,V_shaft(:,y),~,NC_c2] = solver_compressX3_pistoni_explicito_ts(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,beta_3,N_aly(y),mappa_c1,mappa_c2,CINC,T_oil,T_water,m_oil,m_water);
    
    [f] = compressor_dynLIMIT_NW (NC_c2,beta_2,mappa_c2,mc_c2test(y));  %% limiti del compressore 2 %%
    flag2 = -2*isnan(f);
    
    P_shaft(y) = sum(V_shaft(:,y));
    
    if flag2 ~= -2  &&  V_shaft(1,y) > 0  &&  V_shaft(2,y) > 0  &&  P_shaft(y) <= P_disp
    
        eta_best(y) = (V_shaft(1,y)*eta_b1(y) + V_shaft(2,y)*eta_b2(y) + V_shaft(3,y)*eta_b3(y))/P_shaft(y);

    else

        eta_best(y) = 0;
    
    end

end

[eb,ib] = max(eta_best);

if eb == 0

    %soluzione non trovata, dovrò mandare l'energia in rete
    pass  = 0; 
    NC_c1 = 0;
    % mc_c1 = 0;
    % mc_c2 = 0;
    % eta_1 = 0;
    % eta_2 = 0;
    % massa = 0;
    % NC_c2 = 0;
    % Matrix_th     = 0;
    % V_power_shaft = 0;
else

    pass  = 1;
    NC_c1 = N_aly(ib);

    % [mc_c1,mc_c2,eta_1,eta_2,Matrix_th(1:2,:),V_power_shaft(1:2),massa,NC_c2] = solver_compressX2_explicito(f_Eta_c1,f_Nc_c2,f_Eta_c2,beta_1,beta_2,N_aly(ib),mappa_c1,mappa_c2,CINC);

    % graf_e = 1;
    % [mc_t1(i),mc_t2(i),mc_t3(i),eta_1(i),eta_2(i),eta_3(i),V_power_shaft(:,i),massa(i),Nc_T2,Nc_T3,Tout_t(i,:)] = solver_TrenoX3_explicito(f_Eta_c_t1,f_Nc_t2,f_Eta_c_t2,f_Nc_t3,f_Eta_c_t3,beta_1,beta_2,beta_3,N_aly(ib),mappa_c1,Mappa_t2_corr,Mappa_t3_corr,TINC);    
    % P_ob = P_ob0;
    % flag = 0;
    % clear beta_test mc_test eta_b1 eta_b2 eta_b3 V_shaft P_shaft N_aly eta_best

end

if graf_bc3 == 1

    beta_test1(1:Nsol) = beta_1;
    beta_test2(1:Nsol) = beta_2;
    figure(655)
    plot(mappa_c1(:,:,1),mappa_c1(:,:,2))
    xlabel('m_c [kg/s]'); ylabel('\beta []');
    title('Mappa Compressore Centrifugo 1')
    hold on; grid on
    % plot(M_ch_su_line_c1_corr(:,:,2),choke_line_t1_corr(:,2),'k')
    % plot(surge_line_t1_corr(:,1),surge_line_t1_corr(:,2),'r')
    plot(mc_c1test(1:end),beta_test1(1:end),'xk','LineWidth',1)
    plot(mc_c1test(ib),beta_test1(ib),'xr','LineWidth',1)
    hold off
    
    figure(656)
    plot(mappa_c2(:,:,1),mappa_c2(:,:,2))
    xlabel('m_c [kg/s]'); ylabel('\beta []');
    title('Mappa Compressore Centrifugo 2')
    hold on; grid on
    % plot(choke_line_t2_corr(:,1),choke_line_t2_corr(:,2),'k')
    % plot(surge_line_t2_corr(:,1),surge_line_t2_corr(:,2),'r')
    plot(mc_c2test(1:end),beta_test2(1:end),'xk','LineWidth',1)
    plot(mc_c2test(ib),beta_test2(ib),'xr','LineWidth',1)
    hold off               

end

end



