function u0 = initial(x)
    % Parameter definitions
    b = 1.5;
    z = 0.3;
    delta = 0.005;
    alpha = log(2) / (36 * delta^2);
    beta = 10;
    
    % Define base functions
    G = @(x, alpha, z) exp(-alpha * (x - z).^2);
    F = @(x, beta, b) sqrt(max(1 - beta^2 * (x - b).^2, 0));
    
    % Initialize u0 to zero
    u0 = zeros(size(x));
    
    % Vectorized computation for each segment
    % Segment 1: 0.2*pi <= x <= 0.4*pi
    mask1 = (x >= 0.2*pi) & (x <= 0.4*pi);
    u0(mask1) = (1/6) * (G(x(mask1)/pi, alpha, z-delta) + G(x(mask1)/pi, alpha, z+delta) + 4*G(x(mask1)/pi, alpha, z));

    % Segment 2: 1.8 <= x <= 2.5
    mask2 = (x >= 1.8) & (x <= 2.5);
    u0(mask2) = 1.0;
  
    % mask5 = (x >= 1.875) & (x <= 1.8);
    % u0(mask5) = 40*x(mask5) - 75;
    % mask6 = (x >= 2.5) & (x <= 2.525);
    % u0(mask6) = -40*x(mask6) + 101;

    mask5 = (x >= 1.75) & (x <= 1.8);
    u0(mask5) = 20*x(mask5) - 35;
    mask6 = (x >= 2.5) & (x <= 2.55);
    u0(mask6) = -20*x(mask6) + 51;
    
    %  mask5 = (x >= 1.7) & (x <= 1.8);
    % u0(mask5) = 10*x(mask5) - 17;
    % mask6 = (x >= 2.5) & (x <= 2.6);
    % u0(mask6) = -10*x(mask6) + 26;


    % Segment 3: 3 <= x <= 3.8
    mask3 = (x >= 3.) & (x <= 3.8);
    u0(mask3) = 1.0 - abs(2.5* (x(mask3) - 3.4));

    % Segment 4: 1.4*pi <= x <= 1.6*pi
    mask4 = (x >= 1.4*pi) & (x <= 1.6*pi);
    u0(mask4) = (1/6) * (F(x(mask4)/pi, beta, b-delta) + F(x(mask4)/pi, beta, b+delta) + 4*F(x(mask4)/pi, beta, b));

    % Remaining values are already set to 0 by default

end
