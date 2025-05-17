function plotquiver3(pts, normals)
    quiver3(pts(:,1), pts(:,2), pts(:,3), normals(:,1), normals(:,2), normals(:,3), 'r');
end
