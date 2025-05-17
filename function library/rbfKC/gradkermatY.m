function mat=gradkermatY(X,Y)
% create kernel matrices 
% for two point sets X and Y
% corresponding to the full gradient wrt. the Y variable.
% See gradkermatX for details.
mat=-gradkermatX(X,Y);