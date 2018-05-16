clear all
close all
clc

nr_cores     = [0 1 2 3 4 5 6];

% Timings with a total of 50k particles
timings_50k  = [342.74008 83.05219 20.40140 5.11652 1.32374 0.38502 0.11588];
pressure_50k = [0.202134 0.202634 0.205647 0.207849 0.200321 0.219347 0.194617];

% Timings each with 10k particles
timings_10k_each   = [11.804968 12.334865 12.349498 13.947468 15.357635 13.572056 11.247533];
pressures_10k_each = [0.037142 0.082535 0.165876 0.316338 0.567281 1.090733 1.665412];

h1 = figure(1)

subplot(2,1,1)
bar(nr_cores, timings_50k)

text(0:length(timings_50k)-1, timings_50k, num2str(timings_50k'),'verticalalignment','bottom','horizontalalignment','center')
axis([-1 7 0 400])

xlabel 'Nr. of cores (log2)'
ylabel 'Simulation time (s)'
title  'Nr. of MPI cores vs. simulation time with 50k particles'

subplot(2,1,2)
bar(nr_cores, pressure_50k, 'rx--')
set(gca, 'Xtick', nr_cores)

delta = 0.25;
text(0:length(pressure_50k)-1, pressure_50k, num2str(pressure_50k'),'verticalalignment','bottom','horizontalalignment','center')
axis([-1 7 (1-delta).*min(pressure_50k) (1+delta).*max(pressure_50k)])

xlabel 'Nr. of cores (log2)'
ylabel 'Obtained pressure in the box'
title  'Nr. of MPI cores vs. Average pressure in the box with 50k particles'
%Save this plot
print(h1, 'Speedup_50_particles.png', '-dpng')


h2 = figure(2)
subplot(2,1,1)
bar(nr_cores, timings_10k_each)
axis( [min(nr_cores)-1 max(nr_cores)+1 min(timings_10k_each)-10 max(timings_10k_each)+10] )

text(0:length(timings_10k_each)-1, timings_10k_each, num2str(timings_10k_each'),'verticalalignment','bottom','horizontalalignment','center')
xlabel 'Nr. of cores (log2)'
ylabel 'Simulation time (s)'
title  'Nr. of MPI cores vs. simulation time with 10k particles per core'

subplot(2,1,2)
bar(nr_cores,  pressures_10k_each, 'ro--')
set(gca, 'Xtick', nr_cores)

delta = 0.12;
text(0:length(pressures_10k_each)-1, pressures_10k_each, num2str(pressures_10k_each'),'verticalalignment','bottom','horizontalalignment','center')
axis([-1 7 (1-delta).*min(pressures_10k_each) (1+delta).*max(pressures_10k_each)])

xlabel 'Nr. of cores (log2)'
ylabel 'Obtained pressure in the box'
title  'Nr. of MPI cores vs. Average pressure in the box with 10k particles / core'
% Save this plot
print(h2, 'Pressure_increase.png', '-dpng')