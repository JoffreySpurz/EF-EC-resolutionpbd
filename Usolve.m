%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B : Matrice (nP x nP) de la forme A\Mc
% Fc : Vecteur (nP x 1) du laser
% U : matrice de la temperature dans la plaque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ U ] = Usolve(B, Fc, niter)
%%% Condition initiale %%%
U = sparse(nP, niter); 
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Resolution %%%
for i=2:niter
    U(:,k) = B*(U(:,k-1)+Fc);
end
%%%%%%%%%%%%%%%%%%%%%
end

