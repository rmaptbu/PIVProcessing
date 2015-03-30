clear;clc;close all

load('G:\New\Mats\PIV_0604')

npbs=1:35;
%% continuity of spatially resolved fluxes - terms in haematocrit weighted volumetric flow
F_A_rbc_1=(F(:,1,5,1)-F(:,1,1,1)-F(:,1,3,1))./F(:,1,5,1)*100;
F_A_rbc_2=(F(:,1,3,1)-F(:,1,2,1)-F(:,1,4,1))./F(:,1,3,1)*100;
F_A_rbc_t=(F(:,1,5,1)-F(:,1,1,1)-F(:,1,2,1)-F(:,1,4,1))./F(:,1,5,1)*100;

dF_A_rbc_1=1./F(:,1,5,1).*sqrt((dF(:,1,5,1).*(F(:,1,1,1)+F(:,1,3,1))./F(:,1,5,1)).^2+dF(:,1,1,1).^2+dF(:,1,3,1).^2)*100;
dF_A_rbc_2=1./F(:,1,3,1).*sqrt((dF(:,1,3,1).*(F(:,1,2,1)+F(:,1,4,1))./F(:,1,3,1)).^2+dF(:,1,2,1).^2+dF(:,1,4,1).^2)*100;
dF_A_rbc_t=1./F(:,1,5,1).*sqrt((dF(:,1,5,1).*(F(:,1,1,1)+F(:,1,2,1)+F(:,1,4,1))./F(:,1,5,1)).^2+dF(:,1,1,1).^2+dF(:,1,2,1).^2+dF(:,1,4,1).^2)*100;

F_N_rbc_1=(F(npbs,1,5,2)-F(npbs,1,1,2)-F(npbs,1,3,2))./F(npbs,1,5,2)*100;
F_N_rbc_2=(F(npbs,1,3,2)-F(npbs,1,2,2)-F(npbs,1,4,2))./F(npbs,1,3,2)*100;
F_N_rbc_t=(F(npbs,1,5,2)-F(npbs,1,1,2)-F(npbs,1,2,2)-F(npbs,1,4,2))./F(npbs,1,5,2)*100;

dF_N_rbc_1=1./F(npbs,1,5,2).*sqrt((dF(npbs,1,5,2).*(F(npbs,1,1,2)+F(npbs,1,3,2))./F(npbs,1,5,2)).^2+dF(npbs,1,1,2).^2+dF(npbs,1,3,2).^2)*100;
dF_N_rbc_2=1./F(npbs,1,3,2).*sqrt((dF(npbs,1,3,2).*(F(npbs,1,2,2)+F(npbs,1,4,2))./F(npbs,1,3,2)).^2+dF(npbs,1,2,2).^2+dF(npbs,1,4,2).^2)*100;
dF_N_rbc_t=1./F(npbs,1,5,2).*sqrt((dF(npbs,1,5,2).*(F(npbs,1,1,2)+F(npbs,1,2,2)+F(npbs,1,4,2))./F(npbs,1,5,2)).^2+dF(npbs,1,1,2).^2+dF(npbs,1,2,2).^2+dF(npbs,1,4,2).^2)*100;

[e1_A,~,se1_A]=Wmean(abs(F_A_rbc_1),dF_A_rbc_1,1);
[e2_A,~,se2_A]=Wmean(abs(F_A_rbc_2),dF_A_rbc_2,1);
[et_A,~,set_A]=Wmean(abs(F_A_rbc_t),dF_A_rbc_t,1);

[e1_N,~,se1_N]=Wmean(abs(F_N_rbc_1),dF_N_rbc_1,1);
[e2_N,~,se2_N]=Wmean(abs(F_N_rbc_2),dF_N_rbc_2,1);
[et_N,~,set_N]=Wmean(abs(F_N_rbc_t),dF_N_rbc_t,1);

[e1_A se1_A e2_A se2_A et_A set_A ]
[e1_N se1_N e2_N se2_N et_N set_N ]




