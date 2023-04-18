function [bc_type] = bc_type()

global nx;
global ny;

bc_type = [];

% walls
for j = [1,ny]
    for i = 2:nx-1
        k = (j-1)*nx+i;
        bc_type(end+1,:) = [k,1];
    end
end

% far-field (corners included)
for i = [1,nx]
    for j = 1:ny
        k = (j-1)*nx+i;
        bc_type(end+1,:) = [k,2];
    end
end

end

