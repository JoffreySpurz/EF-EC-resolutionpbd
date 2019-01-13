%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A : Matrice (nP x nP) du membre de gauche
% Mc : Matrice (nP x nP) de masse modifiee
% Fc : Vecteur (nP x 1) du laser
% I : liste des points du bord
% U : matrice de la temperature dans la plaque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ U ] = Usolve(A, Mc, Fc, I, U)
%%% Prise en compte des CL %%%
% Nombre de points du bord
nI = length(I);

% CL de Dirichlet
A(I,:) = 0;
A(I,I) = speye(nI);
Mc(I,:) = 0;
Mc(I,I) = speye(nI);
%Fc(I) = 0;% utile ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Resolution %%%
%U = A\(Mc*(U+Fc));% Donnee par les tuteurs
U = A\(Mc*U+Fc);% Bonne formulation ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

