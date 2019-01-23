%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U : matrice de la temperature dans la plaque
% Psi : matrice de psi(x,t)
% B : Matrice (nP x nP) de la forme A\Mc
% nP : Nombre de sommets du maillage
% niter : Nombre d'iterations pour la resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Adjoint ] = adjoint(U,Psi,B,beta_0,beta_1,lambda_0,nP,niter)
%%% Expression des adjoints intermediares
% Pour G : lambda -> u
Adjoint_lambda = adjoint_lambda(U,Psi,B,nP,niter);
% Pour G : beta -> u
Adjoint_beta = adjoint_beta(U,Psi,B,nP,niter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Expression finale : Adjoint = beta_0*beta_1*Adj_beta + lambda_0*Adj_lambda
Adjoint = beta_0*beta_1*Adjoint_beta + lambda_0*Adjoint_lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

