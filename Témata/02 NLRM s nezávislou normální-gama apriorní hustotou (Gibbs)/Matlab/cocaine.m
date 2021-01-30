clear all
close all
clc

%% Pridani cesty k podpurnym funkcim v adresari Support
% adresar pro vlastni funkce (pokud nejsou ve stejnem adresari jako skript)
 addpath('.\Support')

%% Nacteni dat
%56 pozorovani promennych vztahujicich se k prodeji kokainu v
%severovychodni Kalifornii v obodbi 1984-1991. Jedna se o podmnozinu dat
%vyuzitych ve studii
%Caulkins, J.P. a R. Padman (1993): "Quantity Discounts and Quality Premia
%for Illicit Drugs," Journal of the American Statistical Association, 88,
%748-757

% price ... cena za prodany gram kokainu v dané transakci
% quant ... poèet gramù prodaných v dané transakci
% qual  ... kvalita kokainu vyjádøená jako procento èistoty
% trend ... casova promenna, kde 1984=1 a 1991=8

load data_cocaine.mat

