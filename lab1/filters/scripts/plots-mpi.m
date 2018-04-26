clear all
close all
clc

[a,b,c,d] = textread("../results-thresc-mpi.txt", "%s %s %d %f");
threads=c(1:end);
threstimes=d(1:end);
[a,b,c,d] = textread("../results-blurc-mpi.txt", "%s %s %d %f");
blurtimes=d(1:end);

h1 = figure(1)
bar(log2(threads), 1000.*threstimes)
grid on
ylim([0 160])

title('Thresholding filter (MPI) - execution times')
xlabel("Nr cores used (log2)")
ylabel("Execution time (ms)")

# Save the image
print(h1, "thresc-mpi-results.png", "-dpng")

h2 = figure(2)
bar(log2(threads), 1000.*blurtimes, "facecolor", "r")
grid on

title('Blur filter (MPI) - execution times')
xlabel("Nr cores used (log2)")
ylabel("Execution time (ms)")

# Save the image
print(h2, "blurc-mpi-results.png", "-dpng")
