% SET_DIMENSIONS(HANDLE)
% sets orientation, size, etc.

function set_dimensions(handle)
figure(handle);
clf;
orient('landscape');
set(handle,'papertype', '<custom>')
set(handle,'paperunits','centimeters');
set(handle,'papersize',[19 10])
set(handle,'paperposition', [0,0,[19 10]])
end
