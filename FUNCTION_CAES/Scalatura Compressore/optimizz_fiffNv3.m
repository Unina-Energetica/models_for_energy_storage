function Youtput= optimizz_fiffNv3 (X)

global mc_rated_new m_rated beta_rated beta_rated_new Ti pi Tinew pinew


zita = X(1);
% r_n = X(2) ;
r_n = sqrt(Tinew/Ti)*(1/zita);

Youtput(1) = mc_rated_new -m_rated*r_n*(Ti/Tinew)*(pinew/pi)*zita^3;
% Youtput(1) = r_n - sqrt(Tinew/Ti)*(1/zita);
% Youtput(1) = beta_rated_new - (1 + r_n^2*zita^2*(Ti/Tinew)*(beta_rated^eps - 1))^(1/eps); 
% Youtput(1) = beta_rated_new - (1 + (sqrt(Tinew/Ti)*(1/zita))^2*zita^2*(Ti/Tinew)*(beta_rated^eps - 1))^(1/eps); 

end