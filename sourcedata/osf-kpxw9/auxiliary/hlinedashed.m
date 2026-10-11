function hlinedashed(y_pos,col)
if nargin < 2
    col = 'k';
end
range=axis;
hold on
plot([range(1) range(2)],[y_pos y_pos],':','color',col)
hold off
