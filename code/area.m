function [area] = area(xy)

global nx;
global ny;

area_triangles = [];

% Choose a node that is not on the boundaries
ind_x = floor(nx/2);
ind_y = floor(ny/2);
ind = (ind_y - 1)*nx + ind_x;
A = xy(ind,:);

% 1st case
ind_1 = (ind_y - 1)*nx + ind_x - 1;
B = 0.5*(A + xy(ind_1,:));

ind_2 = ind_y*nx + ind_x;
D = 0.5*(A + xy(ind_2,:));

C = A + (2/3)*(0.5*(xy(ind_1,:) + xy(ind_2,:)) - A);

a = norm(B-A);
b = norm(C-B);
c = norm(C-A);
d = norm(D-C);
e = norm(D-A);

s1 = 0.5*(a+b+c);
s2 = 0.5*(c+d+e);
area_triangles(1) = sqrt(s1*(s1-a)*(s1-b)*(s1-c));
area_triangles(2) = sqrt(s2*(s2-c)*(s2-d)*(s2-e));


% 2nd case
ind_3 = ind_y*nx + ind_x;
B = 0.5*(A + xy(ind_3,:));

ind_4 = ind_y*nx + ind_x + 1;
D = 0.5*(A + xy(ind_4,:));

C = A + (2/3)*(0.5*(xy(ind_3,:) + xy(ind_4,:)) - A);

a = norm(B-A);
b = norm(C-B);
c = norm(C-A);
d = norm(D-C);
e = norm(D-A);

s1 = 0.5*(a+b+c);
s2 = 0.5*(c+d+e);
area_triangles(3) = sqrt(s1*(s1-a)*(s1-b)*(s1-c));
area_triangles(4) = sqrt(s2*(s2-c)*(s2-d)*(s2-e));


% 3rd case
ind_5 = ind_y*nx + ind_x + 1;
B = 0.5*(A + xy(ind_5,:));

ind_6 = (ind_y - 1)*nx + ind_x + 1;
D = 0.5*(A + xy(ind_6,:));

C = A + (2/3)*(0.5*(xy(ind_5,:) + xy(ind_6,:)) - A);

a = norm(B-A);
b = norm(C-B);
c = norm(C-A);
d = norm(D-C);
e = norm(D-A);

s1 = 0.5*(a+b+c);
s2 = 0.5*(c+d+e);
area_triangles(5) = sqrt(s1*(s1-a)*(s1-b)*(s1-c));
area_triangles(6) = sqrt(s2*(s2-c)*(s2-d)*(s2-e));


% Construct the area for each cell
for i = 1:nx
    for j = 1:ny
        if (i == 1 && j == 1) || (i == nx && j == ny)
            area((j-1)*nx + i,1) = area_triangles(3) + area_triangles(4) + area_triangles(5) + area_triangles(6);
        elseif (i == 1 && j == ny) || (i == nx && j == 1)
            area((j-1)*nx + i,1) = area_triangles(1) + area_triangles(2);
        elseif ismember(i, [1,nx]) || ismember(j, [1,ny])
            area((j-1)*nx + i,1) = sum(area_triangles);
        else
            area((j-1)*nx + i,1) = 2*sum(area_triangles);
        end
    end
end

end

