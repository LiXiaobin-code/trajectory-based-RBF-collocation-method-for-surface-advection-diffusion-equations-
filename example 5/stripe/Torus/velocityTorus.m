function v = velocityTorus(x,R,r)
% This function contains the exact velocity for the cases of the surface
% advection on the Torus.

% rho = sqrt(x(:,1).^2+x(:,2).^2);
% theta = 1/3*atan(x(:,2)./x(:,1));
% phi = 1/2*atan(x(:,3)./(rho-R));
% rho1 = R + r*cos(2*phi);
% drho1 = -2*r*sin(2*phi);
% 
% v = [-3*rho1.*sin(3*theta) + drho1.*cos(3*theta),...
%     3*rho1.*cos(3*theta)+drho1.*sin(3*theta),...
%     2*r*cos(2*phi)];


[s,rho,~] = cart2pol(x(:,1),x(:,2),x(:,3));
[t, ~] = cart2pol(rho-R, x(:,3));
T1 = [-sin(s) cos(s) 0*s];
T2 = [-sin(t).*cos(s) -sin(t).*sin(s) cos(t)];

v = 3*(R+r*cos(t)).*T1 + 2*r*T2;