%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pi est un sommet du triangle
% a=[a1 a2 a3] matrice 2*3
% b=[b1 b2 b3] matrice 1*3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ a, b ] = base( P1, P2, P3 )
%¨On définit la matrice M
M = [P1' 1; P2' 1; P3' 1];

% On inverse M
M1 = inv(M);

% On définit les ak et bk
a = M1(1:2 , 1:3);
b = M1(3 , 1:3);
end

