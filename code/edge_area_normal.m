function [edge_area_normal] = edge_area_normal(xy)

global nx;
global ny;

edge_area_normal = [];
nu = [];

% horizontal edges
for j = 1:ny
    for i = 1:nx-1
        k = (j-1)*nx+i;
        P = xy(k,:);

        if j == 1
            A = P + (2/3)*(0.5*(xy(k+nx+1,:) + xy(k+1,:)) - P);
            B = 0.5*(P + xy(k+1,:));
            l1 = norm(B-A);
            nu1 = normal_vector(B-A)/norm(normal_vector(B-A));
            edge_area_normal(end+1,:) = l1*nu1;
            
        elseif j == ny
            B = 0.5*(P + xy(k+1,:));
            C = P + (2/3)*(0.5*(xy(k+1,:) + xy(k-nx,:)) - P);
            l2 = norm(C-B);
            nu2 = normal_vector(C-B)/norm(normal_vector(C-B));
            edge_area_normal(end+1,:) = l2*nu2;

        else
            A = P + (2/3)*(0.5*(xy(k+nx+1,:) + xy(k+1,:)) - P);
            B = 0.5*(P + xy(k+1,:));
            C = P + (2/3)*(0.5*(xy(k+1,:) + xy(k-nx,:)) - P);
    
            l1 = norm(B-A);
            l2 = norm(C-B);
    
            nu1 = normal_vector(B-A)/norm(normal_vector(B-A));
            nu2 = normal_vector(C-B)/norm(normal_vector(C-B));
    
            edge_area_normal(end+1,:) = l1*nu1 + l2*nu2;
        end
    end
end

% vertical edges
for j = 1:ny-1
    for i = 1:nx
        k = (j-1)*nx+i;
        P = xy(k,:);

        if i == 1
            B = 0.5*(P + xy(k+nx,:));
            C = P + (2/3)*(0.5*(xy(k+nx,:) + xy(k+nx+1,:)) - P);
            l2 = norm(C-B);
            nu2 = normal_vector(C-B)/norm(normal_vector(C-B));
            edge_area_normal(end+1,:) = l2*nu2;

        elseif i == nx
            A = P + (2/3)*(0.5*(xy(k-1,:) + xy(k+nx,:)) - P);
            B = 0.5*(P + xy(k+nx,:));
            l1 = norm(B-A);
            nu1 = normal_vector(B-A)/norm(normal_vector(B-A));
            edge_area_normal(end+1,:) = l1*nu1;

        else
            A = P + (2/3)*(0.5*(xy(k-1,:) + xy(k+nx,:)) - P);
            B = 0.5*(P + xy(k+nx,:));
            C = P + (2/3)*(0.5*(xy(k+nx,:) + xy(k+nx+1,:)) - P);
    
            l1 = norm(B-A);
            l2 = norm(C-B);
    
            nu1 = normal_vector(B-A)/norm(normal_vector(B-A));
            nu2 = normal_vector(C-B)/norm(normal_vector(C-B));
    
            edge_area_normal(end+1,:) = l1*nu1 + l2*nu2;
        end
    end
end

% diagonal edges
for j = 1:ny-1
    for i = 1:nx-1
        k = (j-1)*nx+i;
        P = xy(k,:);

        A = P + (2/3)*(0.5*(xy(k+nx,:) + xy(k+nx+1,:)) - P);
        B = 0.5*(P + xy(k+nx+1,:));
        C = P + (2/3)*(0.5*(xy(k+nx+1,:) + xy(k+1,:)) - P);

        l1 = norm(B-A);
        l2 = norm(C-B);

        nu1 = normal_vector(B-A)/norm(normal_vector(B-A));
        nu2 = normal_vector(C-B)/norm(normal_vector(C-B));

        edge_area_normal(end+1,:) = l1*nu1 + l2*nu2;
    end
end

end