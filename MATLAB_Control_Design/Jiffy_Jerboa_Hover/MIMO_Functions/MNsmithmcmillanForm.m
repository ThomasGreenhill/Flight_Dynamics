%% MNsmithmcmillanForm.m
%   Calculates the Smith-Mcmillan form of a _symbolic_ transfer function
%   matrix of any size.
%
%   Inputs:
%       A, B, C, D: State-Space Matrices
%   __or__
%       p: polynomial matrix equal to the transfer function matrix G divided by d
%       d: greatest common multiple of denominators of G
%
%   Outputs:
%       UL: Left unimodular matrix
%       Mp: Smith-McMillan form of G
%       UR: Right unimodular matrix
%
%
%   Function Calls:
%       MNsmithForm(), developed by Kristof Pucejdl
%
%
%   Notes: 
%       UL*G*UR = Mp
%
%
%   History:
%       04.07.2021: Created and debugged, TVG
%
%

function [UL, Mp, UR] = MNsmithmcmillanForm(A, B, C, D)

syms s

if nargin == 4
    CT = C.';
    SI = s*eye(size(A));
    G = CT*(SI - A)^-1*B + D;

    d = 1/simplify(det(SI - A));
    p = simplify(CT*adjoint(SI-A)*B + D*det(SI-A));
elseif nargin == 2
    p = A;
    d = B;
    G = p*d;
end

if numel(A) > 25
        fprintf("Large TFM. Computation of Smith-McMillan Form Might Take a While... \n\n")
end

Ms = MNsmithForm(p);

Mp = simplify(d*Ms);

[M, N] = size(G);

if M > N % Over actuated plant
    pt = p.';
    [ULw,~] = hermiteForm(pt, s);
    G = G.';
    Mp = Mp.';
    URw = (ULw*G)\Mp;
    UR = simplify(ULw.');
    UL = simplify(URw.');    
else
    [UL,~] = hermiteForm(p, s);
    UR = simplify((UL*G)\Mp);
    
end

end
