clear all
close all
clc

[nr_threads,timings,iterations,temperatures] = textread('output.txt', "%d %f %d %f");

h1 = figure(1, 'papersize', [800 480])

bar(timings) % The x axis is already from 1 - 16, so no need to specify it

xlabel 'Nr. of threads'
ylabel 'Execution time (s)'
title 'Execution time vs. Nr. of threads used'
axis( [min(nr_threads)-1 max(nr_threads)+1 0 max(timings)+1] )

% Workaround: round the decimals to 2
for i = 1: length(timings)
  temp = num2str(timings(i));
  temp = temp(1:4);
  timings(i) = str2double(temp);
end

text(1:length(timings), timings, num2str(timings),'verticalalignment','bottom','horizontalalignment','center')

% Save the plot
print(h1, 'exec-times.png', '-dpng')