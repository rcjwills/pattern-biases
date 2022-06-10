function newmap = diverging0(color_cell,middle,m,lims)

%white0map = white0(color_cell,middle,m+2);
%newmap = white0map([1:m/2 m/2+3:end],:);

if nargin < 4
    white0map = white0(color_cell,middle,2*m);
else
    white0map = white0(color_cell,middle,2*m,lims);
end
newmap = white0map([1:2:m m+2:2:end],:);
