%% infittitore mappa %%
function [mappa_t1_enh]=infittisci_mappa_compressore(mappa,est)
% clear all
% close all
% clc

% load('Mappa_t1.mat')
% load('Mappa_t1_corr.mat')

% mappa   = Mappa_t1;

global graf

Mc      = mappa(:,:,1);
Beta    = mappa(:,:,2);
Nc      = mappa(:,:,3);
ETA     = mappa(:,:,4);

% Essendo la matrice di partenza di 6 colonne, per equispaziare la matrice
% finale, questa deve avere un numero di colonne pari 6 più multipli di 5 per le
% spaziature, poiché ho 5 spazi intermedi (quindi le possibilità sono: 6, 11, 16, 21,...)
% Nel caso della turbina ho 9 + 8 spazi intermendi

% est = 21;
% est = 25;

[row,col,~] = size(Mc);

Mc_exp   = zeros(row,est); 
Beta_exp = zeros(row,est); 
Nc_exp   = zeros(row,est); 
ETA_exp  = zeros(row,est);

spaz = (est-col)/(col-1);

Mc_exp(:,1:(spaz+1):end)   = Mc; 
Beta_exp(:,1:(spaz+1):end) = Beta; 
Nc_exp(:,1:(spaz+1):end)   = Nc; 
ETA_exp(:,1:(spaz+1):end)  = ETA;

% Nc_exp(:,2:2:end-2) = 0.5*Nc(:,1:end-1)+0.5*Nc(:,2:end);

pos = 1;       %posizione
conta_pos = 1; %counter posizione, mi serve per saltare le colonne già piene dall'inizio
inter = 1;     % conteggio dell'intermezzo per accedere ai valori necessadi di Nc iniziale
for i = 1:est
    
    if i == pos
        
    else
        
        Nc_exp(:,i) = (conta_pos*(Nc(:,inter+1)-Nc(:,inter))/(spaz+1) + Nc(:,inter));
        
        if conta_pos == spaz
            
            pos = pos + (spaz+1);
            conta_pos = 1;
            inter = inter + 1;
            
        else
            
            conta_pos = conta_pos + 1;
            
        end
        
    end
    
end
v_nc = Nc_exp(1,:);
tester = Mc(:,1);

pos = 1;       %posizione
conta_pos = 1; %counter posizione, mi serve per saltare le colonne già piene dall'inizio
inter = 1;     % conteggio dell'intermezzo per accedere ai valori necessadi di Nc iniziale
for k = 1:est
    
    if k == pos
        
    else
        
        for i=1:length(tester)
            
            Beta_exp(i,k) = interp1([Nc(i,inter); Nc(i,inter+1)],[Beta(i,inter) Beta(i,inter+1)],Nc_exp(i,k),'makima');
            Mc_exp(i,k)   = interp1([Nc(i,inter); Nc(i,inter+1)],[Mc(i,inter) Mc(i,inter+1)],Nc_exp(i,k),'makima');
            ETA_exp(i,k)  = interp1([Nc(i,inter); Nc(i,inter+1)],[ETA(i,inter) ETA(i,inter+1)],Nc_exp(i,k),'makima');
            
        end
        
            if conta_pos == spaz

                pos = pos + (spaz+1);
                conta_pos = 1;
                inter = inter + 1;

            else

                conta_pos = conta_pos + 1;

            end
            
        
        
    end
    
end

if graf == 1
    
    figure(100)

    subplot(1,2,1)
    plot(Mc,Beta)
    xlabel('m_c [kg/s]'); ylabel('\beta []')
    hold on; grid on
    contour(Mc,Beta,ETA) 
    colorbar
    hold off

    subplot(1,2,2)
    plot(Mc_exp,Beta_exp)
    xlabel('m_c [kg/s]'); ylabel('\beta []')
    hold on; grid on
    contour(Mc_exp,Beta_exp,ETA_exp) 
    colorbar
    hold off

    pause
    
end
mappa_t1_enh(:,:,1) = Mc_exp; 
mappa_t1_enh(:,:,2) = Beta_exp; 
mappa_t1_enh(:,:,3) = Nc_exp; 
mappa_t1_enh(:,:,4) = ETA_exp;

% mappa_t1_corr_enh(:,:,1) = Mc_exp; 
% mappa_t1_corr_enh(:,:,2) = Beta_exp; 
% mappa_t1_corr_enh(:,:,3) = Nc_exp; 
% mappa_t1_corr_enh(:,:,4) = ETA_exp;

% save('mappa_t1_enh.mat','mappa_t1_enh')
% save('mappa_t1_corr_enh.mat','mappa_t1_corr_enh')


