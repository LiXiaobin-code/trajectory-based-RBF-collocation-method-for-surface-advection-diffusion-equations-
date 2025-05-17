function vernormal = vnormals(vertices, faces)
    % vertices is the list of vertices, size Nx3
    % faces is the list of faces, size Mx3, indexing into vertices
    % vernormal will be the list of vertex normals, size Nx3

    % Initialize the vertex normals array
    vernormal = zeros(size(vertices));
    
    % Compute face normals
    for i = 1:size(faces, 1)
        % Extract vertices of the face
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);

        % verticesectors along two edges of the triangle
        edge1 = v2 - v1;
        edge2 = v3 - v1;

        % Normal of the triangle
        faceNormal = cross(edge1, edge2);
        
        % Add the normal to each vertex's normal sum
        for j = 1:3
            vernormal(faces(i, j), :) = vernormal(faces(i, j), :) + faceNormal;
        end
    end
    
    % Normalize each vertex's normal

    vernormal = vernormal ./ vecnorm(vernormal, 2, 2);

end
