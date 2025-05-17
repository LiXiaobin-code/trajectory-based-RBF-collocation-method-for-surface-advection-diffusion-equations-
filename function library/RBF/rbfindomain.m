function [rbf,dxrbf,evalrbf,lap] = rbfindomain(type)
    switch lower(type)
%         case 'wendland11' 
%             rbf = @(e,x,xi) max(1-e*abs(x-xi),0).^3.*(3*e*abs(x-xi)+1);
%             dxrbf = @(e,x,xi) -12*(x-xi)*e^2.*max(1-e*abs(x-xi),0).^2;
%             evalrbf = @(e,r) max(1-e*r,0).^3.*(3*e*r+1);
%             %drbf = @(e,r) -12*r*e^2.*max(1-e*r,0).^2;
%         case 'wendland12' 
%             rbf = @(e,x,xi) max(1-e*sqrt((x-xi).^2),0).^5.*(8*e^2*(x-xi).^2+5*e*sqrt((x-xi).^2)+1);
%             dxrbf = @(e,x,xi) -14*e^2*(x-xi).*(4*e*abs(x-xi)+1).*max(1-e*abs(x-xi),0).^4;
%             evalrbf = @(e,r) max(1-e*r,0).^5.*(8*e^2*r.^2+5*e*r+1); 
%             %drbf = @(e,r) -14*e^2*r.*max(1-e*r,0).^4.*(4*e*r+1);
        case 'wendland31' 
            rbf = @(e,x,xi) max(1-e*sqrt((x-xi).^2),0).^4.*(4*e*sqrt((x-xi).^2)+1);
            dxrbf = @(e,x,xi) -20*(x-xi)*e^2.*max(1-e*sqrt((x-xi).^2),0).^3;
            evalrbf = @(e,r) max(1-e*r,0).^4.*(4*e*r+1);
            lap=@(e,r) 20*e^2.*(4*e*r-1).*max(1-e*r,0).^2;
        case 'wendland32' 
            rbf = @(e,x,xi) max(1-e*sqrt((x-xi).^2),0).^6.*(35*e^2*(x-xi).^2+18*e*sqrt((x-xi).^2)+3);
            dxrbf = @(e,x,xi) -56*e^2*(x-xi).*(5*e*abs(x-xi)+1).*max(1-e*abs(x-xi),0).^5;
            evalrbf = @(e,r) max(1-e*r,0).^6.*(35*e^2*r.^2+18*e*r+3);
            lap=@(e,r) 56*e^2.*(35*(e*r).^2-4*e*r-1).*max(1-e*r,0).^4;
        case 'quadraticmatern'
            rbf = @(e,x,xi) exp(-e*abs(x-xi)).*(3+3*e*abs(x-xi)+(e*abs(x-xi)).^2);
            dxrbf = @(e,x,xi) -e^2*(x-xi).*exp(-e*abs(x-xi)).*(e*abs(x-xi)+1);
            evalrbf = @(e,r) exp(-e*r).*(3+3*e*r+(e*r).^2);
            lap=@(e,r)  e^2*exp(-e*r).*(-1-e*r+(e*r).^2);
        case 'mq'
            rbf = @(e,x,xi) sqrt(1+(e.*(x-xi)).^2);
            dxrbf = @(e,x,xi) e.^2.*(x-xi)./sqrt(1+(e.*(x-xi)).^2);
            evalrbf = @(e,r) sqrt(1+(e.*r).^2);
            lap=@(e,r) -e.^4.*r.^2./(sqrt(1+(e.*r).^2)).^3+e.^2./sqrt(1+(e.*r).^2);
        case 'gaussian'
            rbf = @(e,x,xi) exp(-e.^2.*(x-xi).^2);
            dxrbf = @(e,x,xi)  -2.*e.^2.*(x-xi).*exp(-e.^2.*(x-xi).^2);  
            evalrbf = @(e,r) exp(-(e.*r).^2);
            lap=@(e,r) -2.*e.^2.*exp(-e.^2.*r.^2)+4.*e.^4.*r.^2.*exp(-e.^2.*r.^2);
    end

    
end

