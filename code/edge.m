function [E] = edge(nx,ny)

global nx;
global ny;

N = (nx-1)*ny+(ny-1)*nx+(nx-1)*(ny-1);

E = zeros(N,2);

k = 1;

% horizontal edges 
for j = 1:ny 
    for i= 1:nx-1
        E(k,1) = (j-1)*nx+i;
        E(k,2) = (j-1)*nx+i+1;
        k = k+1;
    end
end

%vertical edges
for j = 1:ny-1 
    for i= 1:nx
        E(k,1) = (j-1)*nx+i;
        E(k,2) = (j-1)*nx+i+nx;
        k = k+1;
    end
end

% diagonal edges
for j = 1:ny-1
    for i = 1:nx-1
        E(k,1) = (j-1)*nx+i;
        E(k,2) = (j-1)*nx+i+nx+1;
        k = k+1;
    end
end

end

