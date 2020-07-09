function to_gid(init_coord, curr_coord, corresp, result, ele_type)
%{
This function creates the GID post-processing files

Args:
init_coord: inital nodal coordinates
curr_coord: nodal coordinates of the deformed configuration
corresp:    correspondence matrix
results:    stress array
ele_type:   number of nodes per element(20 = quadratic)

Returns: GID input files for post-processing
%}

n_nod = size(init_coord, 1);    % number of nodes
n_el  = size(corresp,    1);    % number of elements

%% Results file:

fid = fopen('salida.post.res', 'w');
if fid == -1
    error('no se puede abrir el archivo: %s', 'salida.txt')
end

% header:
fprintf(fid, 'GID Post Results File 1.0\n');
fprintf(fid, '\n');
fprintf(fid, 'GaussPoints  "GP_HEXAHEDRA_8"        Elemtype   HEXAHEDRA\n');
fprintf(fid, 'Number of Gauss Points:   8\n');
fprintf(fid, 'Natural Coordinates   : Internal\n');
fprintf(fid, 'End GaussPoints\n');
fprintf(fid, '\n');

% the stress in each Gauss point are written:
fprintf(fid, 'Result "STRESSES" "MECHANICAL"              1.0 Matrix OnGaussPoints "GP_HEXAHEDRA_8"\n');
fprintf(fid, 'ComponentNames "Sxx" "Syy" "Szz" "Sxy" "Syz" Sxz"\n');
fprintf(fid, 'Values \n');
for ee = 1:n_el
    for pt = 1:8
        if pt == 1
            fprintf(fid, '%d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n', ee, ...
                    result{ee, pt}(1), result{ee, pt}(2), result{ee, pt}(3),...
                    result{ee, pt}(4), result{ee, pt}(5), result{ee, pt}(6));
        else
            fprintf(fid, '%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n',...
                    result{ee, pt}(1), result{ee, pt}(2), result{ee, pt}(3),...
                    result{ee, pt}(4), result{ee, pt}(5), result{ee, pt}(6));
        end
    end
    
    if ee == n_el
        fprintf(fid, 'End Values \n');
    end
end

% it's writted the nodal displacements:
fprintf(fid, 'Result "disp" "disp"              1.0 Vector OnNodes "GP_HEXAHEDRA_8"\n');
fprintf(fid, 'ComponentNames "Dx" "Dy" "Dz"\n');
fprintf(fid, 'Values \n');
disp = curr_coord - init_coord;
for nd = 1:n_nod
    fprintf(fid, '%d %15.8f %15.8f %15.8f\n', nd,...
            disp(nd, 1), disp(nd, 2), disp(nd, 3));
    
    if nd == n_nod
        fprintf(fid, 'End Values \n');
    end
end
fclose(fid);

%% Mesh file:

fid = fopen('salida.post.msh', 'w');
if fid == -1
    error('no se puede abrir el archivo: %s', 'salida.post.msh');
end

% header:
if ele_type == 20
    fprintf(fid, 'MESH dimension 3 ElemType HEXAHEDRA Nnode 20\n');
end

% nodal coordinates in the deformed configuration:
fprintf(fid, 'Coordinates\n');
for nd = 1:n_nod
    fprintf(fid, '%d %15.8f %15.8f %15.8f\n', nd,...
            curr_coord(nd, 1), curr_coord(nd, 2), curr_coord(nd, 3));    
    
    if nd == n_nod
        fprintf(fid, 'End Coordinates \n');
    end
end

% correspondence matrix according to the element type:
fprintf(fid, 'Elements\n');
if ele_type == 20
    for el = 1:n_el
        fprintf(fid, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',...
                el,  corresp(el, 1),  corresp(el, 3),  corresp(el, 5),...
                     corresp(el, 7), corresp(el, 13), corresp(el, 15),...
                    corresp(el, 17), corresp(el, 19),  corresp(el, 2),...
                     corresp(el, 4),  corresp(el, 6),  corresp(el, 8),...
                     corresp(el, 9), corresp(el, 10), corresp(el, 11),...
                    corresp(el, 12), corresp(el, 14), corresp(el, 16),...
                    corresp(el, 18), corresp(el, 20), 1);

        if el == n_el
            fprintf(fid, 'End Elements \n');
        end
    end
end
fclose(fid);
end