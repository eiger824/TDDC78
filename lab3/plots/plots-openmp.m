clear all
close all
clc

[nr_threads,timings,iterations,temperatures] = textread('output.txt', "%d %f %d %f");

h1 = figure(1)

bar(timings) % The x axis is already from 1 - 16, so no need to specify it

xlabel 'Nr. of threads'
ylabel 'Execution time (s)'
title 'Execution time vs. Nr. of threads used'

text(0:length(timings)-1, timings, num2str(timings'),'verticalalignment','bottom','horizontalalignment','center')

% Save the plot
print(h1, 'exec-times.png', '-dpng')