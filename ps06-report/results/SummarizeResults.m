% SummarizeResults produces graphs from measured data
% outputs both on screen (gnuplot) and into pdfs

% study case, OpenMP or MPI
set_dimensions(1);
omp = load('omp.txt');
mpi = load('mpi.txt');
ompRes = average_runs(omp, 2, 5);
mpiRes = average_runs(mpi, 1, 5);
ompRes(:,2) = time_to_speedup(ompRes(:,2));
mpiRes(:,2) = time_to_speedup(mpiRes(:,2));
x = 1:max([size(ompRes),size(mpiRes)]);
p = plot(ompRes(:,1), ompRes(:,2), 'o-', mpiRes(:,1), mpiRes(:,2), 'o-', x, x, '-');
set(p, 'linewidth', 3);
axis tight
legend('OpenMP', 'MPI', 'reference', 'Location', 'Northwest');
title('OpenMP / MPI speedup');
xlabel('number of nodes / threads');
ylabel('measured speedup');
print -landscape -dpdf studycase.pdf

% convergence rate
set_dimensions(2);
conv = load('convergence.txt');
convRes = average_runs(conv, 3, 4);
x = 0:4;
reference = (10 .^ x)';
reference(:,2) = 10 .^ (-x*2);
p = loglog(convRes(:,1), convRes(:,2), 'o-', reference(:,1), reference(:,2), '-');
set(p, 'linewidth', 3);
legend('measured error', 'reference', 'Location', 'Northeast');
title('Convergence rate test');
xlabel('problem size');
ylabel('maximal error');
print -dpdf convergence.pdf

% OpenMP/MPI tradeoff
set_dimensions(3);
c = load('c.txt');
cRes = average_runs(c, 1, 5);
serialRun = 794.7076;
cRes(:,2) = serialRun ./ cRes(:,2);
p = plot(cRes(:,1), cRes(:,2), 'o-');
set(p, 'linewidth', 3);
title('OpenMP / MPI tradeoff on 36 processors');
xlabel('number of MPI nodes (out of 36)');
ylabel('speedup');
print -landscape -dpdf tradeoff.pdf

% different problem sizes
set_dimensions(4);
hold on;
b = load('b.txt');
b(:,1) = b(:,1) .* b(:,2);
b(:,2) = [];
color = ['b','r','g'];
i = 1;
for ps = unique(b(:,2))',
  ids = find(b(:,2) == ps);
  bRes = average_runs(b(ids,:), 1, 4);
  bRes(:,2) = time_to_speedup(bRes(:,2));
  p = plot(bRes(:,1), bRes(:,2), strcat(color(1,i),'o-'));
  set(p, 'linewidth', 3);
  i = i + 1;
end 
p = plot(bRes(:,1), bRes(:,1), 'c-');
set(p, 'linewidth', 3);
axis tight
legend('8192', '16384', '32768', 'reference', 'Location', 'Northwest');
title('Speedup in different problem sizes');
xlabel('number of processors');
ylabel('speedup');
hold off;
print -dpdf problems.pdf
