clear all
close all
clc

[a,b,c,d] = textread("results-thresc-mpi.txt", "%d %d %d %f");
threads=c(2:end);
threstimes=d(2:end);
[a,b,c,d] = textread("results-blurc-mpi.txt", "%d %d %d %f");
blurtimes=d(2:end);

h1 = figure(1)
bar(log2(threads), 1000.*threstimes)
grid on

title('Thresholding filter (MPI) - execution times')
xlabel("Nr cores used")
ylabel("Execution time (ms)")

# Save the image
print(h1, "thresc-mpi-results.png", "-dpng")

h2 = figure(2)
bar(threads, 1000.*blurtimes, "facecolor", "r")
grid on

title('Blur filter (MPI) - execution times')
xlabel("Nr cores used")
ylabel("Execution time (ms)")

# Save the image
print(h2, "blurc-mpi-results.png", "-dpng")
