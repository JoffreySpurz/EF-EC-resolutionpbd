clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Load %%%%%%%%%%
% Avec (1) ou sans (0) fissure ?
isfissure = 1;

% Resolution
if isfissure == 1
    load('Probleme_direct_avec_fissure_rectangulaire','U', 'mesh1', 'niter_laser')
else
    load('Probleme_direct_sans_fissure','U', 'mesh1', 'niter_laser')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Affichage %%%%%%%%%%%
niter = length(U(1,:));
pas = 1;

%%% Champ de temperature
% Avec laser
for i = 1:pas:niter_laser
    subplot(2,2,1)
    
    %%% Affichage
    mesh1.surf(U(:,i))
    %caxis([0 4]*10^7)
    pause(0.1)
end

% % Sans laser
% for i = niter_laser+1:pas:niter
%     subplot(2,2,2)
%     
%     %%% Affichage
%     mesh1.surf(U(:,i))
%     %caxis([0 4]*10^7)
%     pause(0.1)
% end
% 
% %%% Gradient de temperature
% % Avec laser
% for i = 1:pas:niter_laser
%     subplot(2,2,3)
%     
%     % Calcul de la norme du gradient
%     gu = mesh1.grad(U(:,i));
%     ngu = mesh1.norm(gu,2); 
%     
%     %%% Affichage
%     mesh1.surf(ngu)
%     %caxis([0 4]*10^7)
%     pause(0.1)
% end
% 
% % Sans laser
% for i = niter_laser+1:pas:niter
%     subplot(2,2,4)
%     
%     % Calcul de la norme du gradient
%     gu = mesh1.grad(U(:,i));
%     ngu = mesh1.norm(gu,2); 
%     
%     %%% Affichage
%     mesh1.surf(ngu)
%     %caxis([0 4]*10^7)
%     pause(0.1)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%