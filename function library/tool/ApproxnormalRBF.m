function [ Normal ] = ApproxnormalRBF(Z,num_pt,idx_Global_Local )
global RBFscale
[nZ,dim]=size(Z);
Normal=zeros(nZ,3);
switch idx_Global_Local
    %% by Global pts
    case 'G'
        
         Mat = kermat( Z, Z );
        lam = Mat\ ones( nZ, 1) ;
        % lam'*Mat*lam
        
        tmp=gradkermatX(Z, Z,RBFscale);
        size( tmp)
        for ii=1:3

            Global_normal(:,ii)= tmp(:,:,ii)*lam;
        end
        Global_normal=Global_normal/norm(Global_normal);
        Normal=Global_normal;
        %% RBF normal from local stencil
    case 'L'
        
          Mat = kermat( Z, Z );
        lam = Mat\ ones( nZ, 1) ;
        % lam'*Mat*lam
        tmp=gradkermatX(Z, Z,RBFscale);
        for ii=1:3
            Global_normal(:,ii)= tmp(:,:,ii)*lam;
        end
        Global_normal=Global_normal/norm(Global_normal);
        
        loc_index  = knnsearch(Z ,Z ,'k',num_pt);
        for i=1:length(Z) 
        Z_loc = Z(loc_index(i,:)',:);
        Mat =kermat( Z_loc, Z_loc );
        lam = Mat\ ones( length(loc_index(i,:)), 1);
        % lam'*Mat*lam
        tmp=gradkermatX(Z(i,:), Z_loc,RBFscale);
        for ii=1:3
            Local_normal(ii,1)= tmp(1,:,ii)*lam;
        end
        Normal(i,:) = Local_normal/norm(Local_normal)*sign( Local_normal(1)*Global_normal(i,1));
        end
      
end
end


