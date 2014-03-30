% OUT = TIME_TO_SPEEDUP(DATA)
% transforms timing data to speedup data
% DATA is vector
% speedup computed as a ratio to highest mean (assumed to be 1 proc setup)
% largest data value is considered serial time

function [out] = time_to_speedup(data)
serial = max(data);
out = serial ./ data;
end
