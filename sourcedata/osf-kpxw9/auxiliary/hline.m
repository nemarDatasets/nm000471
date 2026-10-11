function hline(y_pos,col)
% vline(x_pos,col)
if nargin < 2
    col = 'k';
end
range=axis;
hold on
plot([range(1) range(2)],[y_pos y_pos],col)
hold off
