function [v_reshaped] = reshape_array(v)

global nx;
global ny;

for j = 1:ny
    for i = 1:nx
        v_reshaped(j,i) = v((j-1)*nx+i);
    end
end

end

