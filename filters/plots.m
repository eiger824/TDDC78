clear all
close all
clc

[a,b,c,d] = textread("results-thresc-pthreads.txt");
threads=c(2:end);
threstimes=d(2:end);
[a,b,c,d] = textread("results-blurc-pthreads.txt");
blurtimes=d(2:end);

h1 = figure(1)
bar(log2(threads), 1000.*threstimes)
grid on

title('Thresholding filter (PThreads) - execution times')
xlabel("Nr threads used (log2)")
ylabel("Execution time (ms)")

# Save the image
print(h1, "thresc-pthreads-results.png", "-dpng")

h2 = figure(2)
bar(log2(threads), 1000.*blurtimes, 'r')
grid on

title('Blur filter (PThreads) - execution times')
xlabel("Nr threads used (log2)")
ylabel("Execution time (ms)")

# Save the image
print(h2, "blurc-pthreads-results.png", "-dpng")