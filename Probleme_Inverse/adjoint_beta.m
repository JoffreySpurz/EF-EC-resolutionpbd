%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U : matrice de la temperature dans la plaque
% Psi : matrice de psi(x,t)
% B : Matrice (nP x nP) de la forme A\Mc
% nP : Nombre de sommets du maillage
% niter : Nombre d'iterations pour la resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Adjoint ] = adjoint_beta(U,Psi,B,nP,niter)
%%% Obtention de Phi(x,T-t)
% Resolution de l'equation de la chaleur
Phi = Usolve(B, Psi, nP, niter);
% Changement de variables t->T-t
Phi = fliplr(Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calcul de l'adjoint de la différentielle : Adjoint = Dt(u) Phi(x,T-t)
% Troncature sur U(t) pour obtenir U(t-dt)
Utronc = [zeros(nP,1) U(:,1:niter-1)];
% adjoint de la differentielle
Adjoint = (U-Utronc)/dt .* Phi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

