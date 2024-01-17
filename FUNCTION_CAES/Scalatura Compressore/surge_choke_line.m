%% SOUGE LINE %% 
% creazione dati relative alla surge e choke line
function [M_ch_su_line_t1]=surge_choke_line(mappa)
% close all
% clear all
% clc 

% load('Mappa_t1')
% load('Mappa_t1_corr')

% mappa = Mappa_t1;

global graf

Mc   = mappa(:,:,1);
Beta = mappa(:,:,2);
Nc   = mappa(:,:,3);

v_NC = Nc(1,:);

NC_input = linspace(v_NC(1),v_NC(end),1500);

for i=1:length(NC_input)
    
    nc = NC_input(i);
    S_v_NC = v_NC>nc;
    L_v_NC = v_NC<nc;

    NC_sup = S_v_NC.*v_NC;
    NC_inf = L_v_NC.*v_NC;

    NC_sup(NC_sup==0) = NaN;
    [x_sup,y_sup] = min(NC_sup);
    [x_inf,y_inf] = max(NC_inf);

    %% surge section %%

    beta_ext1 = Beta(1,y_inf);
    beta_ext2 = Beta(1,y_sup);
    mc_ext1   = Mc(1,y_inf);
    mc_ext2   = Mc(1,y_sup);

    V_BETA_MAX(i) = interp1([v_NC(y_inf) v_NC(y_sup)],[beta_ext1 beta_ext2 ],nc);
    M_betaM(i)    = interp1([beta_ext1 beta_ext2 ],[mc_ext1  mc_ext2],V_BETA_MAX(i));


    %% choke section %%

    beta_ext1 = Beta(end,y_inf);
    beta_ext2 = Beta(end,y_sup);
    mc_ext1   = Mc(end,y_inf);
    mc_ext2   = Mc(end,y_sup);

    V_BETA_min(i) = interp1([v_NC(y_inf) v_NC(y_sup)],[beta_ext1 beta_ext2 ],nc);
    M_betam(i)  = interp1([beta_ext1 beta_ext2 ],[mc_ext1  mc_ext2], V_BETA_min(i));


end

if graf == 1
    
    figure(999)
    plot(Mc(:,1:end),Beta(:,1:end))
    xlabel('m_c [kg/s]'); ylabel('\beta []')
    title('Mappa Corretta Secondo Compressore')
    hold on; grid on
    plot(M_betam,V_BETA_min,'k')
    plot(M_betaM,V_BETA_MAX,'r')
    legend(num2str(v_NC(1)),num2str(v_NC(2)),num2str(v_NC(3)),num2str(v_NC(4)),num2str(v_NC(5)),num2str(v_NC(6)),num2str(v_NC(7)),num2str(v_NC(8)),num2str(v_NC(9)),num2str(v_NC(10)),num2str(v_NC(11)),num2str(v_NC(12)),num2str(v_NC(13)),num2str(v_NC(14)),num2str(v_NC(15)),num2str(v_NC(16)),num2str(v_NC(17)),num2str(v_NC(18)),num2str(v_NC(19)),num2str(v_NC(20)),num2str(v_NC(21)),'choke','surge')
    contour(Mc,Beta,mappa(:,:,4))
    colorbar
    hold off

    pause
    
end

M_ch_su_line_t1(:,:,1) = [M_betaM',V_BETA_MAX',NC_input']; %surge line
M_ch_su_line_t1(:,:,2) = [M_betam',V_BETA_min',NC_input']; %choke line

% M_ch_su_line_t1_corr(:,:,1) = [M_betaM',V_BETA_MAX',NC_input']; %surge line
% M_ch_su_line_t1_corr(:,:,2) = [M_betam',V_BETA_min',NC_input']; %choke line

% save('M_ch_su_line_t1.mat','M_ch_su_line_t1')
% save('M_ch_su_line_t1_corr.mat','M_ch_su_line_t1_corr')




