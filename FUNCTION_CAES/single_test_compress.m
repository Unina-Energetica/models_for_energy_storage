function [Mc_fit,beta_fit] = single_test_compress(nc,mappa)
Mc      = mappa(:,:,1);
Beta    = mappa(:,:,2);
Nc      = mappa(:,:,3);
% PShaft  = mappa(:,:,4);
% ETAa    = mappa(:,:,5);

v_NC = Nc(1,:);
S_v_NC = v_NC>nc;
L_v_NC = v_NC<nc;

% if nc > mappa(1,end,3)
%     nc = mappa(1,end,3)
% elseif nc < mappa(1,1,3)
%     nc = mappa(1,1,3)
% end 

NC_sup = S_v_NC.*v_NC;
NC_inf = L_v_NC.*v_NC;

[x_inf,y_inf] = max(NC_inf);
NC_sup(NC_sup==0) = NaN;
[x_sup,y_sup] = min(NC_sup);

beta_s = Beta(:,y_sup);
beta_i = Beta(:,y_inf);
Mc_s   = Mc(:,y_sup);
Mc_i   = Mc(:,y_inf);

for i=1:length(beta_s)

    appo_beta(i) = interp1([Nc(i,y_inf); Nc(i,y_sup)],[beta_i(i) beta_s(i)],nc,'makima');
    appo_Mc(i)   = interp1([Nc(i,y_inf); Nc(i,y_sup)],[Mc_i(i) Mc_s(i)],nc,'makima');

end

Mc_fit = appo_Mc;
beta_fit = appo_beta;

