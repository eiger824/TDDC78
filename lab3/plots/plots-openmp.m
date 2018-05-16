clear all
close all
clc

nr_threads = [1 2 4 8 16];
timings    = [4.1334 2.2808 1.3298 0.8115 0.8240];

h1 = figure(1)

bar(log2(nr_threads), timings)

xlabel 'Nr. of threads (log2)'
ylabel 'Execution time (s)'
title 'Execution time vs. Nr. of threads used'

text(0:length(timings)-1, timings, num2str(timings'),'verticalalignment','bottom','horizontalalignment','center')

% Save the plot
print(h1, 'exec-times.png', '-dpng')