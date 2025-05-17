clear
close all


%uu=load('redbell_stripe_5048_t10_tf5281');

matX=load('redcell_6208');
X=matX.p;
 FaceX=matX.t;
 
 
  figure
trisurf(FaceX,X(:,1),X(:,2),X(:,3),X(:,3))%uu.
% shading interp

axis equal
zlim([-1,1])
axis off
colormap(jet)
[cmax,cmin]=clim;
% cmax=-.24;cmin=.23;
climit = get(gca, 'CLim');
    set(gca,'CLim',[cmax,cmin]);
    %set(gca,'Visible','off')
    set(gca,'color','none')
    h.EdgeColor='none';
    h.LineStyle='none';
%     h.FaceColor=[1,0,0];
      shading interp ;
lighting phong;
light
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 15;
h.BackFaceLighting = 'unlit';
set(gca,'FontSize',10)