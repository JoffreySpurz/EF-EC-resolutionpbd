clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Load %%%%%%%%%%
% Resolution
load('Probleme_direct_sans_fissure','U', 'mesh1')
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Affichage %%%%%%%%%%%
niter = length(U(1,:));
pas = 1;

%%% Avec laser
for i = 1:pas:niter
    %%% Affichage
    mesh1.surf(U(:,i))
    %caxis([0 0.05])
    pause(0.01)
end

% % On indique que le laser va etre coupe
% mesh1.surf(zeros(nP,1));
% pause(1);

% %%% Sans laser
% for i = 11:niter
%     %%% Affichage
%     mesh1.surf(U)
%     %caxis([0 0.05])
%     pause(0.1)
%     
%     %%% Iteration
%     U = Usolve(A, M, zeros(nP,1), I, U);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%