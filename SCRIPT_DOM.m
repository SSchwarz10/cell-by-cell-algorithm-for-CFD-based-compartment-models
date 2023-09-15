%% SCRIPT for determining the degree of mixing
clear all;
clc;

% Liu (2012) Age distribution and the degree of mixing in continuous flow
% stirred tank reactors
load('DATA_RAW.mat','M1','M2','M3','M4','M5','Volume','tau_CFD','M1_grad','outlet_massflow');
n_cells = length(M1);
c0_M1_M2_M3_M4_M5 = [[1:n_cells].' M1 M2 M3 M4 M5];
% delete all cells, which are not connected to the outlet (Fluent and
% Matlab enumeration!)
c0_M1_M2_M3_M4_M5 = c0_M1_M2_M3_M4_M5(outlet_massflow(:,1)+1,:);

% Flow Average of the moments at the outlet
M1_e = sum(c0_M1_M2_M3_M4_M5(:,2).*outlet_massflow(:,3))/sum(outlet_massflow(:,3));
M2_e = sum(c0_M1_M2_M3_M4_M5(:,3).*outlet_massflow(:,3))/sum(outlet_massflow(:,3));
M3_e = sum(c0_M1_M2_M3_M4_M5(:,4).*outlet_massflow(:,3))/sum(outlet_massflow(:,3));
M4_e = sum(c0_M1_M2_M3_M4_M5(:,5).*outlet_massflow(:,3))/sum(outlet_massflow(:,3));
M5_e = sum(c0_M1_M2_M3_M4_M5(:,6).*outlet_massflow(:,3))/sum(outlet_massflow(:,3));

M1_sqr_e = sum(c0_M1_M2_M3_M4_M5(:,2).^2.*outlet_massflow(:,3))/sum(outlet_massflow(:,3));

% variance, CoV and skewness at the end of the reactor
variance_e = M2_e - M1_e^2;
CoV_e = sqrt(variance_e)/M1_e;
central_M3_e = M3_e - 3*M2_e*M1_e^1 + 2*M1_e^3;
skewness_e = (central_M3_e/M1_e^3)^(1/3);

% Volume average of M1
M1_Vavg = sum(M1.*Volume)/sum(Volume);
% Volume average of M2
M2_Vavg = sum(M2.*Volume)/sum(Volume);
% Volume average of M1^2
M1_sqr_Vavg = sum(M1.^2.*Volume)/sum(Volume);

% volume averaged numerical diffusivity (Liu, 2012, Age distribution in the
% Kenics Static Micromixer with convection and Diffusion)
M1_grad_magn = sqrt(sum(M1_grad.^2,2));
M1_grad_sqr_Vavg = sum(M1_grad_magn.^2.*Volume)/sum(Volume);
D_eff = (2*tau_CFD*M1_Vavg - M1_sqr_e)/(2*tau_CFD*M1_grad_sqr_Vavg);

% CoV regarding age, mean age and degree of mixing Zwietering and
% Danckwerts (defined as ratio of vaiances w.r.t to volume averaged age)
CoV_age_ZD = sqrt(M2_Vavg - M1_Vavg^2)/M1_Vavg;
CoV_mean_age_ZD = sqrt(M1_sqr_Vavg - M1_Vavg^2)/M1_Vavg;
J_ZD = (M1_sqr_Vavg - M1_Vavg^2)/(M2_Vavg - M1_Vavg^2);

% CoV regarding age, mean age and degree of mixing Liu (defined as ratio of
% vaiances w.r.t to mean residence time)

CoV_age_L = sqrt(M2_Vavg - 2*M1_Vavg*tau_CFD + tau_CFD^2)/tau_CFD;
CoV_mean_age_L = sqrt(M1_sqr_Vavg - 2*M1_Vavg*tau_CFD + tau_CFD^2)/tau_CFD;
J_L = (M1_sqr_Vavg - 2*M1_Vavg*tau_CFD + tau_CFD^2) / (M2_Vavg - 2*M1_Vavg*tau_CFD + tau_CFD^2);

% frequeny function of the CFD and ideal reactor models
[NumBins,BinEdges,BinIndex] = histcounts(M1,'BinMethod','fd');
BinMiddle = (BinEdges(2:end) + BinEdges(1:end-1))/2;
BinWidth = BinEdges(2) - BinEdges(1);   % all bins with same width
M1_V = [M1 Volume];
volume_fraction = zeros(length(NumBins),1);
for i = 1:length(NumBins)
    volume_fraction(i) = sum(M1_V(BinIndex==i,2))/sum(M1_V(:,2));
end

g_CFD = volume_fraction/BinWidth;
g_CFD_scaled = g_CFD*tau_CFD;

g_PFR_scaled = ones(length(NumBins),1);
g_PFR_scaled(BinMiddle/tau_CFD>1) = 0;

g_LTR_scaled = tau_CFD^2./(BinMiddle.').^2/4;
g_LTR_scaled(BinMiddle < tau_CFD/2) = 1;

g_CSTR_scaled = [1 0;
          1 2];
plot(BinMiddle/tau_CFD,g_CFD_scaled,'LineWidth',1.5);
hold on
grid on
plot(BinMiddle/tau_CFD,g_PFR_scaled,'LineWidth',1.5);
plot(BinMiddle/tau_CFD,g_LTR_scaled,'LineWidth',1.5);
plot(g_CSTR_scaled(:,1),g_CSTR_scaled(:,2),'LineWidth',1.5);
xlabel('M_{1,\theta} in -')
ylabel('g(M_{1,\theta}) in -')
xlim([0 BinEdges(end)/tau_CFD])
legend('CFD','PFR','LTR','CSTR')
xlim([0 2])
ylim([0 1.2])
title('Frequency function g(M_{1,\theta}) for CFD and different CPT models')