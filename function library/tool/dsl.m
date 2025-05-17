function dsl(m,x1,y1,x2)

% drawslop: (no input) draw a line and print slope in log-log plot based on two mouse clicks
% drawslop(m): draw a line with slope m starting from first mouse click and end at the
% x-locoation of the second 
%
% (updated on 11/20/2023) to allow 4 inputs to plot without ginput

hold on
if nargin==0
    g=ginput(2);
    x1=g(1,1); y1=g(1,2); x2=g(2,1); y2= g(2,2); %y1*(x2/x1)^(m);
    m = log(y2/y1)/log(x2/x1);
    plot( [x1,x2], [y1,y2],'k--')
    text( 0.5*(x1+x2),y2, sprintf('slope %0.1f',m),'FontSize',16)
%     text( 0.1*(x1+x2) + max(x1,x2), y2, sprintf('slope %0.1g',m),'FontSize',16)
elseif nargin==1
    g=ginput(2);
    x1=g(1,1); y1=g(1,2); x2=g(2,1); y2= y1*(x2/x1)^(m);
    plot( [x1,x2], [y1,y2],'k--')
    text( 0.5*(x1+x2),y2, sprintf('slope %0.1d',m),'FontSize',16)
elseif nargin==4
    y2= y1*(x2/x1)^(m);
    plot( [x1,x2], [y1,y2],'k--', 'HandleVisibility', 'off')
    text( 0.4*(x1+x2),y2/1e3, sprintf('slope %0.1d',m),'FontSize',16)
else
    warning('Incorrect inputs')
end
