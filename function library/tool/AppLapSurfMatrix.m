function [Lap] = AppLapSurfMatrix(X,Surfcetype,num_pt)
%only for cyclide by KPM
global RBFtype;% RBFtype='ms';  % see frbf.m

global RBFscale

% out normal vector
if num_pt>0
    %%PCA
%[ Normal ] = ApproxnormalPCA( X, num_pt );

%%Global 
idx_Global_Local='G';
[ Normal ] =ApproxnormalRBF(X,num_pt,idx_Global_Local );

%%Local
% idx_Global_Local='L';
% [ Normal ] =ApproxnormalRBF(X,num_pt,idx_Global_Local );

else
    % analytic normal
func_normal=cell(3,1);

syms x y z
switch Surfcetype
    
    case 1 %unit sphere
        func_normal{1}=@(x,y,z)x;
        func_normal{2}=@(x,y,z)y;
        func_normal{3}=@(x,y,z)z;
    case 3 %Torus
        
        func_normal{1}=matlabFunction(1/((4*(x^2+y^2+z^2+8/9)*x-8*x)^2+(4*(x^2+y^2+z^2+8/9)*y-8*y)^2+16*(x^2+y^2+z^2+8/9)^2*z^2)^(1/2)*(4*(x^2+y^2+z^2+8/9)*x-8*x));
        func_normal{2}=matlabFunction(1/((4*(x^2+y^2+z^2+8/9)*x-8*x)^2+(4*(x^2+y^2+z^2+8/9)*y-8*y)^2+16*(x^2+y^2+z^2+8/9)^2*z^2)^(1/2)*(4*(x^2+y^2+z^2+8/9)*y-8*y));
        func_normal{3}=matlabFunction(4/((4*(x^2+y^2+z^2+8/9)*x-8*x)^2+(4*(x^2+y^2+z^2+8/9)*y-8*y)^2+16*(x^2+y^2+z^2+8/9)^2*z^2)^(1/2)*(x^2+y^2+z^2+8/9)*z);
end

Nx=func_normal{1}(X(:,1),X(:,2),X(:,3));
Ny=func_normal{2}(X(:,1),X(:,2),X(:,3));
Nz=func_normal{3}(X(:,1),X(:,2),X(:,3));
R=sqrt(Nx.^2+Ny.^2+Nz.^2);
Nx=Nx./R;
Ny=Ny./R;
Nz=Nz./R;
Normal=[Nx(:),Ny(:),Nz(:)];
end
Pij=cell(3,3);

gradmat=gradkermatX(X,X,RBFscale);
G=zeros(length(X),length(X),3);

   Lap=zeros(length(X),length(X));
    for i1=1:3
        for j=1:3
            if (i1==j)
           Pij{j,i1}=diag(1-Normal(:,i1).*Normal(:,j));
            else
                Pij{j,i1}=diag(-Normal(:,i1).*Normal(:,j));
            end
        end
        G(:,:,i1)=(Pij{1,i1}*gradmat(:,:,1)+Pij{2,i1}*gradmat(:,:,2)+Pij{3,i1}*gradmat(:,:,3));% 
        Lap=Lap+G(:,:,i1)/matrixgen( X, X,[0,0,0,1], RBFscale)*G(:,:,i1);
    end
 

end

