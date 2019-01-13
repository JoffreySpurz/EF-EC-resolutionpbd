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
niter = 15; % Nombre d'iterations pour la resolution
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
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv *10^(-12); % Puissance volumique du laser utilisee
% Second membre
Fc = dt*S*mesh1.P1('x.^2+y.^2<0.2^2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%

%%% Avec laser
for i = 1:10
    %%% Affichage
    mesh1.surf(U)
    caxis([0 0.05])
    pause(0.1)
    
    %%% Iteration
    U = Usolve(A, M, Fc, I, U);
end

% On indique que le laser va etre coupe
mesh1.surf(zeros(nP,1));
pause(1);

%%% Sans laser
for i = 11:niter
    %%% Affichage
    mesh1.surf(U)
    caxis([0 0.05])
    pause(0.1)
    
    %%% Iteration
    U = Usolve(A, M, zeros(nP,1), I, U);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Affichage au temps final %%%%%%%%%
% mesh1.surf(U)
% caxis([0 1000])
% 
% figure;
% mesh1.plot

% figure;
% mesh1.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%