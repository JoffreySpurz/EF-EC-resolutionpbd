%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U : matrice de la temperature dans la plaque
% Psi : matrice de psi(x,t)
% B : Matrice (nP x nP) de la forme A\Mc
% nP : Nombre de sommets du maillage
% niter : Nombre d'iterations pour la resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Adjoint ] = adjoint_lambda(U,Psi,B,nP,niter)
%%% Obtention de Phi(x,T-t)
% Resolution de l'equation de la chaleur
Phi = Usolve(B, Psi, nP, niter);
% Changement de variables t->T-t
Phi = fliplr(Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calcul de l'adjoint de la differentielle : Adjoint = - grad(u).grad( Phi(x,T-t) )
% Gradient de u(x,t) : Matrice (np x 2niter) Les niter premieres colonnes sont Gradx puis on a Grady
GradU = mesh1.grad(U);
% Gradient de phi(x,T-t)
GradPhi = mesh1.grad(Phi);
% adjoint de la differentielle
Adjoint = -(GradU(:,1:niter).*GradPhi(:,1:niter) + GradU(:,niter+1:2*niter).*GradPhi(:,niter+1:2*niter));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

