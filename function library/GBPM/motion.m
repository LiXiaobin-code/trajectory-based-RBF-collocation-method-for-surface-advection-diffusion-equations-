function ctrs = motion(ctrs,V,n,T,dt)
% This function moves the closest points according to a motion law
% using the velocity in the normal and tangential components:
% 
% v = V*n + T

ctrs = ctrs + dt*(bsxfun(@times,V,n) + T);