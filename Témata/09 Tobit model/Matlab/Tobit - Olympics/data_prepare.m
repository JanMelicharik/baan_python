% Priprava dat

clear
close all
clc

data_raw = readtable('olympics.xlsx');

data_olympics = table2struct(data_raw);

save olympics.mat data_olympics