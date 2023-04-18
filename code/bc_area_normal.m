function [bc_area_normal] = bc_area_normal(bc_type)

global nx;
global ny;
global delta_x;
global delta_y;

bc_area_normal = zeros(length(bc_type),2);

for j = 1:ny
    for i = 1:nx
        k = (j-1)*nx+i;

        if ismember(k, bc_type(:,1))
            index = find(k == bc_type(:,1));
    
            % walls
            if bc_type(index,2) == 1
                % walls (not corners)
                if j == 1 && ~ismember(k,[1,nx,nx*(ny-1)+1,ny*nx])
                    nu = delta_x*[0,-1];
                elseif j == ny && ~ismember(k,[1,nx,nx*(ny-1)+1,ny*nx])
                    nu = delta_x*[0,1];
                end
    
            % far-field
            elseif bc_type(index,2) == 2
                if i == 1 && ~ismember(k,[1,nx,nx*(ny-1)+1,ny*nx])
                    nu = [-1,0]*delta_y;
                elseif i == nx && ~ismember(k,[1,nx,nx*(ny-1)+1,ny*nx])
                    nu = [1,0]*delta_y;
                % corners
                elseif i == 1 && j == 1
                    nu = 0.5*(delta_x*[0,-1] + delta_y*[-1,0]);
                elseif i == 1 && j == ny
                    nu = 0.5*(delta_x*[0,1] + delta_y*[-1,0]);
                elseif i == nx && j == 1
                    nu = 0.5*(delta_x*[0,-1] + delta_y*[1,0]);
                elseif i == nx && j == ny
                    nu = 0.5*(delta_x*[0,1] + delta_y*[1,0]);
                end
            end

            bc_area_normal(index,:) = nu;

        end
    end
end


end

