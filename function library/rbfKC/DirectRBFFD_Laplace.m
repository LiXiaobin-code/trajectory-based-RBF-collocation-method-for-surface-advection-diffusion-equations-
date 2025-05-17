function lambda=DirectRBFFD_Laplace(ctrs_h,ctrs,poly)


[N,d]=size(ctrs);



kmat=kermat(ctrs,ctrs);
b=laplacekermat(ctrs_h,ctrs)';


pmat=[ctrs, ones(N,1)];


A=[kmat,pmat;pmat', zeros(d+1,d+1)];
B=[b;zeros(d+1,1)];

switch poly
case 1
lambda=A\B;
lambda=lambda(1:N)';
case 0
lambda=kmat\b;
lambda=lambda';
end


