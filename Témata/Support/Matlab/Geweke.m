function results = Geweke(theta)
% Gewekova konvergenèní diagnostika
% vyuzivajici funkci momentg.m z
% LeSageho ekonometrickeho toolboxu
% theta ... rozmer S1 x k
% (S pocet vzorku, k pocet parametru)
    
    [S1,~] = size(theta);
        
    smpl_A = round(0.1*S1); %prvnich 10% vzorku
    smpl_C = round(0.6*S1)+1; %poslednich 40% vzorku
    
    %NSE pro vzorek A
    pom = momentg(theta(1:smpl_A,:));
    mean_A = [pom.pmean]';
    nse_A = [pom.nse1]';
    
    %NSE pro vzorek B
    pom = momentg(theta(smpl_C+1:S1,:));
    mean_C = [pom.pmean]';
    nse_C = [pom.nse1]';
    
    %Gewekova konvergenèní diagnostika
    CD = (mean_A-mean_C)./(nse_A+nse_C);
    
    %ulozeni vysledku
    results.CD = CD;
    
end

