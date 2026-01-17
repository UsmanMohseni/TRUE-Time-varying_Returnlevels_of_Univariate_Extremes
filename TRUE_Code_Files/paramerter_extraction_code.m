% This code extracts the parameters for median, 10 percentile and 90
% percentile.

clc
clear
close all
%% Load the file
% Please add the "non_EMD_20_evaluation_new.mat" from the Result_2025
% folder for all the stations in 1_Narmada.
% The example shows for hoshangabad station.
% for copying the address, click on the file and press ctrl+shift+C
load("E:\Mayank\1_Narmada\3_hoshangabad\Result_2025\non_EMD_20_evaluations_new.mat");
% x= parameters for specific covariate model
% put x=1 for Annual Precipitation
% put x=2 for LTA
% put x=7 for AP+LTA

x=1;
%%
% All variable names are self representative
% PF=Peak Flow
% PD=Peak Duration
% PV=Peak Volume
% pctle= percentile

% Note : if any PARA_"anyname" is of size 1x14 , it means that flood
% characterstic is stationary, hence parameter have to be calculated
% directly using param function. 
A=median(squeeze(PARA_FINAL_PF1(:,:,x)));
B=prctile(squeeze(PARA_FINAL_PF1(:,:,x)),10);
C=prctile(squeeze(PARA_FINAL_PF1(:,:,x)),90);
PARA_PF_median = A(~isnan(A(:)));
PARA_PF_10prctle = B(~isnan(B(:)));
PARA_PF_90prctle = C(~isnan(C(:)));


D=median(squeeze(PARA_FINAL_PD1(:,:,x)));
E=prctile(squeeze(PARA_FINAL_PD1(:,:,x)),10);
F=prctile(squeeze(PARA_FINAL_PD1(:,:,x)),90);
PARA_PD_median = D(~isnan(D(:)));
PARA_PD_10prctle = E(~isnan(E(:)));
PARA_PD_90prctle = F(~isnan(F(:)));

G=median(squeeze(PARA_FINAL_PV1(:,:,x)));
H=prctile(squeeze(PARA_FINAL_PV1(:,:,x)),10);
I=prctile(squeeze(PARA_FINAL_PV1(:,:,x)),90);
PARA_PV_median = G(~isnan(G(:)));
PARA_PV_10prctle = H(~isnan(H(:)));
PARA_PV_90prctle = I(~isnan(I(:)));