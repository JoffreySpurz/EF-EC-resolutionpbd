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
mesh1 = Mesh(domain1,dx,'save');

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

% Plaque
rho_plaque = 7850;% masse volumique (Kg/m^3) : valeur de l'acier
cp_plaque = 444;% capacite thermique massique (J/K/Kg) : valeur de l'acier
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
lambda_plaque = 50.2; % conductivite thermique (W/m/K) : valeur de l'acier 

% Iteration
niter = 15; % Nombre d'iterations pour la resolution
dt =10^(-6); % Pas en temps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Matrices d'iterations %%%%%%%%
%%% Raideur
% Param�tre : Rho*cp sur chaque triangle
LAMBDA = mesh1.P0(1,2);
% Matrice de raideur
Kc = dt*mesh1.stiffness(LAMBDA);%matriceK(P, T, Lambda, dt);


%%% Masse
% Param�tre : Conductivite thermique (lambda) sur les triangles
BETA = mesh1.P0(1);
% Matrice de masse
M = mesh1.mass(BETA);
%mesh1.mass(Beta.*ones(nT,1));


%%% Membre de gauche
% Probleme AX(n) = M( X(n-1) + F )
A = M + Kc;

%%% Prise en compte des CL %%%
% Nombre de points du bord
nI = length(I);

% CL de Dirichlet
A(I,:) = 0;
A(I,I) = speye(nI);
Mc(I,:) = 0;
Mc(I,I) = speye(nI);

% Probleme X(n) = A\M ( X(n-1) + F ) = B ( X(n-1) + F )
B = A\M;


%%% Laser
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv *10^(-5); % Puissance volumique du laser utilisee
% Second membre
Fc = dt*S*mesh1.P1('exp(-(x.^2+y.^2)/0.1)');% Gaussien
%Fc=dt*S*mesh1.P1('x.^2+y.^2<0.2^2');% Cercle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%
%%% Avec laser
U = Usolve(B, Fc, nP, niter);

%%% Save
save('Probleme_direct_sans_fissure','U', 'mesh1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%