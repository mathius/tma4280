% OUT = AVERAGE_RUNS(IN, COLX, COLY)
% in table format: nodes, threads, problem size, error, time elapsed
% out table format: columns defined by COLX, mean of COLY in corresponding group
% WARNING: COLX can be only 1-3 to work properly!

function [out] = average_runs(in, colx, coly)
in(:,6) = in(:,3)*10000 + in(:,2)*100 + in(:,1);
ids(:,2) = unique(in(:,6));
ids(:,1) = 1:size(ids);
for id = ids(:,1)',
	num = ids(id,2);
	in(find(in(:,6) == num), 7) = id;
end
out = unique(in(:,colx));
out(:,2) = accumarray(in(:,7), in(:, coly), [], @mean);
end
