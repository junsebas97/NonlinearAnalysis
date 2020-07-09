function to_gid(init_coord, curr_coord, corresp, result, ele_type)
%{
This function creates the GID post-processing files

Args:
init_coord: inital nodal coordinates
curr_coord: nodal coordinates of the deformed configuration
corresp:    correspondence matrix
results:    stress array
ele_type:   number of nodes per element(4 = bi-linear FE, 8 = quadratic FE)

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
fprintf(fid, 'GaussPoints  "GP_Quadrilateral"        Elemtype   Quadrilateral \n');
fprintf(fid, 'Number of Gauss Points:   4\n');
fprintf(fid, 'Natural Coordinates   : Internal\n');
fprintf(fid, 'End GaussPoints\n');
fprintf(fid, '\n');

% the stress in each Gauss point are written:
fprintf(fid, 'Result "STRESSES" "MECHANICAL"              1.0 Matrix OnGaussPoints "GP_Quadrilateral"\n');
fprintf(fid, 'ComponentNames "Sxx" "Syy" "Sxy"\n');
fprintf(fid, 'Values \n');
for ee = 1:n_el
    for pt = 1:4
        if pt == 1
            fprintf(fid, '%d %15.8f %15.8f %15.8f\n', ee, ...
                    result{ee, pt}(1), result{ee, pt}(2), result{ee, pt}(3));
        else
            fprintf(fid, '%15.8f %15.8f %15.8f\n',...
                    result{ee, pt}(1), result{ee, pt}(2), result{ee, pt}(3));
        end
    end
    
    if ee == n_el
        fprintf(fid, 'End Values \n');
    end
end

% it's writted the nodal displacements:
fprintf(fid, 'Result "disp" "disp"              1.0 Vector OnNodes "GP_Quadrilateral"\n');
fprintf(fid, 'ComponentNames "Dx" "Dy" \n');
fprintf(fid, 'Values \n');
disp = curr_coord - init_coord;
for nd = 1:n_nod
    fprintf(fid, '%d %15.8f %15.8f\n', nd, disp(nd, 1), disp(nd, 2));
    
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
if ele_type == 4
    fprintf(fid, 'MESH dimension 3 ElemType Quadrilateral Nnode 4\n');
elseif ele_type == 8
    fprintf(fid, 'MESH dimension 3 ElemType Quadrilateral Nnode 8\n');
end

% nodal coordinates in the deformed configuration:
fprintf(fid, 'Coordinates\n');
for nd = 1:n_nod
    fprintf(fid, '%d %15.8f %15.8f %15.8f\n', nd,...
            curr_coord(nd, 1), curr_coord(nd, 2), 0.0);    
    
    if nd == n_nod
        fprintf(fid, 'End Coordinates \n');
    end
end

% correspondence matrix according to the element type:
fprintf(fid, 'Elements\n');
if ele_type == 4
    for el = 1:n_el
        fprintf(fid, '%d %d %d %d %d %d\n', el, corresp(el, 1),...
                corresp(el, 2), corresp(el, 3), corresp(el, 4), 1);

        if el == n_el
            fprintf(fid, 'End Elements \n');
        end
    end
    
elseif ele_type == 8
    for el = 1:n_el
        fprintf(fid, '%d %d %d %d %d %d %d %d %d %d\n', el, corresp(el, 1), ...
                corresp(el, 3), corresp(el, 5), corresp(el, 7), ...
                corresp(el, 2), corresp(el, 4), corresp(el, 6), ...
                corresp(el, 8), 1);

        if el == n_nod
            fprintf(fid, 'End Elements \n');
        end
    end
end
fclose(fid);
end