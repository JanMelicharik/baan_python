function boot_index = my_boot(y,S,w)
%Vazeny bootstrap
% index ... S*1 vektor bootstrapovych indexu
% y ... N x 1 vektor bootstrapovaneho vzorku
% S ... velikost bootstrapoveho vyberu
% w ... bootstrapove vahy (optional)

    N = length(y);

    if nargin < 3
     w = 1/N;
    end

    %vektor indexu
    boot_index = zeros(S,1);

    %vytvoreni bodu intervalu, kde sirka dilcich intervalu reprezentuje
    %stanovenou vahu
    w_int = cumsum(w);

    %nastaveni "seedu" generatoru nahodnych cisel
    %rng(sum(clock)*sum(date));
    
    for s=1:S
        %rand*w_int(end) ...nahodne cislo z uniformniho rozdeleni na intervalu
        %(0;w_int(end)) a ziskani indexu jako suma intervalu splnujici
        %podminku ze hodnota krajni meze intervalu je mensi nez toto nahodne
        %vygenerovane cislo
        boot_index(s) = sum(w_int<rand*w_int(end))+1;
    end


end

