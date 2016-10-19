clear all 
close all
NX = 60;
NY = 220;

load results_tau_0_day.out
sw_no_tau = reshape( results_tau_0_day(:,2),NX,NY);
clear results_tau_0_day

load results_tau_15_day.out
sw_yes_tau = reshape( results_tau_15_day(:,2),NX,NY);
clear results_tau_15_day

figure
subplot(1,2,1)
imagesc(sw_no_tau)
title('Sw at 5 years: tau=0')
subplot(1,2,2)
imagesc(sw_yes_tau)
title('Sw at 5 years: tau=15 day')

figure
imagesc(sw_no_tau - sw_yes_tau)
title('Difference in Sw (tau=0 or 15 day)')
colorbar

