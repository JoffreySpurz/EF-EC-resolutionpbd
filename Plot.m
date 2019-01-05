clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Maillage %%%%%%%%%%
% 1 Unit square
domain1 = Domain('square');

% mesh
dx = 0.5; % taille d'un element
mesh1 = Mesh(domain1,dx);

% Liste des sommets
P = mesh1.nodes;

% Liste des triangles
T = mesh1.triangles;

% Liste des points du bords
I = mesh1.boundary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Parametres %%%%%%%%%%
% Nombre de sommets
nP = length(P); 

% Nombre de triangles
nT = length(T);

% Condition initiale
U = sparse(nP, 1); 

% Plaque
Beta = beta( nT); % Vecteur Rho*cp sur chaque triangle
Lambda = lambda(nT); % Vecteur Conductivite thermique (lambda) sur les triangles
Flaser = flaser(nT); % Vecteur laser sur les triangles

% Iteration
niter = 20; % Nombre d'iterations pour la resolution
dt = 0.1; % Pas en temps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Matrices d'iterations %%%%%%%%
%%% Raideur
% Matrice de raideur
Kc = matriceK(P, T, Lambda, dt);


%%% Masse
% Matrice de masse
Mc = mesh1.mass(Beta.*ones(nT,1));


%%% Membre de gauche
% Probleme AX(n) = MX(n-1) + F
A = Mc + Kc;


%%% Laser
% Second membre
Fc = vecteurF( P, T, Flaser, dt );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%
for i = 1:niter
    %%% Affichage
    mesh1.surf(U)
    caxis([0 1000])
    pause(0.1)
    
    %%% Iteration
    U = Usolve(A, Mc, Fc, I, U);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Affichage au temps final %%%%%%%%%
mesh1.surf(U)
caxis([0 1000])

figure;
mesh1.plot

% figure;
% mesh1.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%