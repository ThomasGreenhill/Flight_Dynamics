function G = sym2tf(g)
%{
SYM2TF Symbolic Transfer Function to Numerical Transfer Function

Conversion

Syntax:  G = sym2tf(g)

Inputs:
   g - Symbolic Transfer Function representation

Outputs:
   G - Numeric Transfer Function representation

Revision History:
   06.09.2020: Created and Debugged, TVG

Credit:
   This program is based off Oskar Vivero Osornio's work but updated and
   improved for R2020a

Example: 
   syms p
   g11=(p + 2)/(p^2 + 2*p + 1);
   g12=(p - 1)/(p^2 + 5*p + 6);
   g21=(p - 1)/(p^2 + 3*p + 2);
   g22=(p + 2)/(p + 1);
   g=[g11 g12; g21 g22];
   G=sym2tf(g)
%}

% g = simplify(g);
[n,m]=size(g);
for i=1:n
    for j=1:m
        [num,den]=numden(g(i,j));
        num_n=sym2poly(num);
        den_n=sym2poly(den);
        G(i,j)=tf(num_n,den_n);
    end
end
