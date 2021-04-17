function g = tf2sym(G)
%
% TF2SYM Numerical Transfer Function to Symbolic Transfer Function**.
% conversion.
% 
% Syntax:  g = tf2sym(G)
% 
% Inputs:
%    G - Numerical Transfer Function
% 
% Outputs:
%    g - Symbolic Transfer Function
% 
% Revision History:
%    06.09.2020: Created and Debugged, TVG
% 
% Credit:
%    This program is based off Oskar Vivero Osornio's work but updated and
%    improved for R2020a
% 
% Example: 
%    g11=tf([1 2],[1 2 1]);
%    g12=tf([1 -1],[1 5 6]);
%    g21=tf([1 -1],[1 3 2]);
%    g22=tf([1 2],[1 1]);
%    G=[g11 g12; g21 g22];
%    g=tf2sym(G);
%    pretty(g)
%

p=sym('s');
[n,m]=size(G);
g=zeros(n,m);
g=sym(g);
for i=1:n
    for j=1:m
        [num,den]=tfdata(G(i,j),'v');
        l=length(den);
        order=l;
        for k=1:l
            A(order,1)=p^(k-1);
            order=order-1;
        end
        n=num*A;
        d=den*A;
        g(i,j)=n/d;
        clear A;
    end
end

g=simplify(g);
