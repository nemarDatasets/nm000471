function plot_contour(x_axis_, y_axis_, cont_matrix, linespec,linewidth, corner_spec)


if nargin < 4
   linespec = 'k-';
   linewidth = 1;
end
if nargin < 6
   corner_spec = 1;
end

if numel(x_axis_) ~= size(cont_matrix,2) || numel(y_axis_) ~= size(cont_matrix,1) 
   error('dimenson mismatch between x or y axis and contour matrix') 
end

[x,y] = meshgrid(x_axis_, y_axis_);
switch corner_spec
    case 1
        x = interp2(x, 4);
        y = interp2(y, 4);
        contourlines = interp2(cont_matrix, 4);
        hold on
        contour(x,y,contourlines,1,linespec,'LineWidth',linewidth);
    case 2
        x = interp2(x, 2);
        y = interp2(y, 2);
        contourlines = interp2(cont_matrix, 2, 'nearest');
        dx = mean(diff(x(1, :)));
        dy = mean(diff(y(:, 1)));
        hold on
        contour(x+dx/2,y+dy/2,contourlines,1,linespec,'LineWidth',linewidth);
end