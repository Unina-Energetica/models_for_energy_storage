function [media, mediana, gaussiana, deviazione_standard] = Statistica(vettore)
    % Calcolo della media
    media = mean(vettore);
    
    % Calcolo della mediana
    mediana = median(vettore);
    
    % Calcolo della distribuzione gaussiana
    pd = fitdist(vettore, 'Normal');
    gaussiana = struct('Mean', pd.mu, 'StandardDeviation', pd.sigma);
    
    % Calcolo della deviazione standard
    deviazione_standard = std(vettore);
end