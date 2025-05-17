function X = pickpointscircle(n, varargin)

    if isempty(varargin)
        radius = 1;
    else
        radius = varargin{1}; 
    end
    

    theta = linspace(0, 2*pi, n+1)'; 
    theta(end) = []; 

    X = [cos(theta) * radius, sin(theta) * radius];

end