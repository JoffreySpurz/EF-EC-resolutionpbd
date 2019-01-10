clear; close all;

tic

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
dx = 0.1; % taille d'un element
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

% Iteration
niter = 1000; % Nombre d'iterations pour la resolution
dt = 0.1; % Pas en temps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Matrices d'iterations %%%%%%%%
%%% Raideur
% Paramètre : Rho*cp sur chaque triangle
LAMBDA = mesh1.P0(1,2);
% Matrice de raideur
Kc = dt*mesh1.stiffness(LAMBDA);%matriceK(P, T, Lambda, dt);


%%% Masse
% Paramètre : Conductivite thermique (lambda) sur les triangles
BETA = mesh1.P0(1);
% Matrice de masse
M = mesh1.mass(BETA);
%mesh1.mass(Beta.*ones(nT,1));


%%% Membre de gauche
% Probleme AX(n) = MX(n-1) + F
A = M + Kc;


%%% Laser
% Second membre
Fc = dt*mesh1.P1('x.^2+y.^2<0.2^2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%

%%% Avec laser
for i = 1:10
    %%% Iteration
    U = Usolve(A, M, Fc, I, U);
end

%%% Sans laser
for i = 11:niter
    %%% Iteration
    U = Usolve(A, M, zeros(nP,1), I, U);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc