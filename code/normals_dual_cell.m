function [normals] = normals_dual_cell(i,j, edge_area_normal, bc_area_normal, bc_type)

global nx;
global ny;

normals = zeros(7,2);
N_1 = ny*(nx-1);
N_2 = N_1 + nx*(ny-1);

if i <= nx-1
    %right
    normals(3,:) = edge_area_normal((j-1)*nx+i-(j-1),:);
end
if i >= 2
    %left
    normals(6,:) = -edge_area_normal((j-1)*nx+i-1-(j-1),:);
end

if j <= ny-1
    %up
    normals(1,:) = edge_area_normal(N_1 + (j-1)*nx+i,:);
end
if j >= 2
    %down
    normals(4,:) = -edge_area_normal(N_1 + (j-1)*nx+i-nx,:);
end

if i <= nx-1 && j <= ny-1
    %up-right
    normals(2,:) = edge_area_normal(N_2 + (j-1)*nx+i-(j-1),:);
end
if i >= 2 && j >= 2
    %down-left
    normals(5,:) = -edge_area_normal(N_2 + (j-1)*nx+i-(j-1)-(nx-1)-1,:);
end

if j == 1 || j == ny || i == 1 || i == nx
    k = (j-1)*nx+i;
    index = find(k == bc_type(:,1));
    normals(7,:) = bc_area_normal(index,:);
end


end

