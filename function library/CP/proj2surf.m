function proj_pts = proj2surf(pts, surftype)
    % Input parameters:
    % pts - Nx3 matrix, where each row represents the coordinates of a point [x, y, z]
    % surftype - string specifying the type of surface ('Sphere', 'Torus', 'Bretzel2', etc.)

    % Output parameters:
    % proj_pts - Nx3 matrix, containing the coordinates of each point projected onto the surface

    % Number of points
    n = size(pts, 1);

    % Initialize output matrix
    proj_pts = zeros(n, 3);

    % Set optimization options
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    % Loop over each point
    for i = 1:n
        % Current point
        pt = pts(i, :);

        % Initial guess
        init_guess = pt;

        % Call fmincon
        [proj_pt, ~] = fmincon(@(x) dist_fun(x, pt), ...
                               init_guess, [], [], [], [], [], [], ...
                               @(x) surf_constr(x, surftype), opts);

        % Store projected point
        proj_pts(i, :) = proj_pt;
    end
end

function dist = dist_fun(x, pt)
    % Calculate distance function
    dist = sum((x - pt).^2);
end

function [c, ceq] = surf_constr(x, surftype)
    % Surface equation constraint
    c = [];  % Inequality constraint (not needed)
    ceq = surf_eq(x, surftype);  % Equality constraint, describes the surface
end

function val = surf_eq(p, surftype)
    % Define surface equation based on the specified surface type
    switch surftype
        case 'Sphere'
            % Equation for a sphere centered at (0,0,0) with radius 1
            val = sum(p.^2) - 1;

        case 'Torus'
            % Equation for a torus with major radius R and minor radius r
            R = 1;
            r = 1/3;
            val = (sum(p.^2) + R^2 - r^2)^2 - 4*R^2*(p(1)^2 + p(2)^2);

        case 'Bretzel2'
            % Equation for the 'Bretzel2' surface
            val = (-p(1)^4 + p(1)^2 - p(2)^2)^2 + p(3)^2/2 - (1 + sum(p.^2))/40;

        case 'Cyclide'
            % Equation for a cyclide surface
            val = (sum(p.^2) - 1 + 19^2/100)^2 - 4*(2*p(1) + sqrt(4 - 19^2/100))^2 - (4*(19*p(2))^2)/100;

        case 'CPD surface'
            % Equation for a complex polynomial surface
            val = sqrt((p(1)-1)^2 + p(2)^2 + p(3)^2) * sqrt((p(1)+1)^2 + p(2)^2 + p(3)^2) * ...
                  sqrt(p(1)^2 + (p(2)-1)^2 + p(3)^2) * sqrt(p(1)^2 + (p(2)+1)^2 + p(3)^2) - 1.1;

        case 'orthocircle'
            % Equation for an orthocircle surface
            val = ((p(1)^2 + p(2)^2 - 1)^2 + p(3)^2) * ...
                  ((p(2)^2 + p(3)^2 - 1)^2 + p(1)^2) * ...
                  ((p(1)^2 + p(3)^2 - 1)^2 + p(2)^2) - 0.075^2 * (1 + 3 * (p(1)^2 + p(2)^2 + p(3)^2));

        case 'Peanut'
            % Equation for a peanut-shaped surface
            val = ((2*p(1)-1)^2 + 4*p(2)^2 + 4*p(3)^2) * ...
                  ((2*p(1)+1)^2 + 4*p(2)^2 + 4*p(3)^2) - 1.2;

        otherwise
            error('Unsupported surface type');
    end
end
