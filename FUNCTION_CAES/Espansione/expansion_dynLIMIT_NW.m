function [f] = expansion_dynLIMIT_NW (nc,beta,mappa,mass_in)
% nota la mappa di fuzionamento del compressore
% si entra con il numero di giri e il rapporto di compressione
% derinisce la portata corretta a quell'istante

if nc < mappa(1,1,3) | nc > mappa(1,end,3)
    
    f = NaN;
    
else
    
    [Mc_fit,beta_fit] = single_test_compress(nc,mappa);

    %scelgo così i valori perché ai limiti le curve tendono a essere
    %orizzontali per le portate minime e verticali per le portate massime
    beta_max_ist = beta_fit(end);
    beta_min_ist = min(beta_fit);
    Mc_max_ist = max(Mc_fit);
    Mc_min_ist = Mc_fit(1);  
    
%     figure(999)
%     plot(mappa(:,:,1),mappa(:,:,2))   
%     hold on   
%     grid on
%     plot(Mc_fit,beta_fit,'b--o')
%     plot(mass_in,beta,'k*')
%     xlabel('m_c [kg/s]'); ylabel('\beta [ ]')
%     title('Mappa Scalata Corretta')
%     hold off 
    
    if    mass_in<Mc_min_ist | beta<beta_min_ist | beta>beta_max_ist  | mass_in>Mc_max_ist % punto di funzionamento innammissiile
%        Mc_ist = 0;
        f = NaN;         
        % beta_max_ist
        % beta_min_ist
        % beta
        % Mc_max_ist
        % Mc_min_ist
        % mass_in

    else

%         Mc_ist = mass_in;
        f = 1;

    end

end

end







