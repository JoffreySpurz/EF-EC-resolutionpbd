%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U : matrice de la temperature dans la plaque
% U : matrice de psi(x,t)
% B : Matrice (nP x nP) de la forme A\Mc
% nP : Nombre de sommets du maillage
% niter : Nombre d'iterations pour la resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ grad ] = adjoint_beta(U,Psi,B,nP,niter)
%%% Obtention de Phi(x,T-t)
% Resolution de l'equation de la chaleur
Phi = Usolve(B, Psi, nP, niter);
% Changement de variables t->T-t
Phi = fliplr(Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calcul du gradient : grad = Dt(u) Phi(x,T-t)
grad = zeros(nP,niter);
for k = 2:niter
    grad(:,k)= (U(:,k)-U(:,k-1))/dt .* Phi(:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%
end

