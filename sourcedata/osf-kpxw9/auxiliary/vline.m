function vline(x_pos,col)
% vline(x_pos,col)
if nargin < 2
    col = 'k';
end
range=axis;
% ymin=range(3);
% ymax=range(4);
hold on
plot([x_pos x_pos],[range(3) range(4)],col)
hold off
