function [xy] = xy(x_min, x_max, y_min, y_max)

global nx;
global ny;
global x_min;
global x_max;
global y_min;
global y_max;
global delta_x;
global delta_y;

xy = [];

for j = 1:ny
    for i = 1:nx
        xy(end+1,:) = [x_min + (i-1)*delta_x, y_min + (j-1)*delta_y];
    end
end

end