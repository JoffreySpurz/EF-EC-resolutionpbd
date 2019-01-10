clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Maillage %%%%%%%%%%
% Fissure carree
pos_x = -0.5;
pos_y = -0.75;
epaisseur = 0.25;
longueur = 1;

% 1 Unit square
nodes = [-1 -1; 1 -1 ; 1 1; -1 1;pos_x pos_y ; pos_x+longueur pos_y; pos_x+longueur pos_y+epaisseur ;pos_x pos_y+epaisseur];
edges = {{1,2} {2,3} {3,4} {4,1} {5,6,'lr'} {6,7,'lr'} {7,8,'lr'} {8,5,'lr'}}; 
domain1 = Domain(nodes,edges);

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

% Plaque
Beta = beta( nT); % Vecteur Rho*cp sur chaque triangle
Lambda = lambda(nT); % Vecteur Conductivite thermique (lambda) sur les triangles

% Iteration
niter = 15; % Nombre d'iterations pour la resolution
dt = 0.1; % Pas en temps
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
% Probleme AX(n) = MX(n-1) + F
A = M + Kc;


%%% Laser
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv * pi*(0.2/10^5)^2; % Puissance du laser applique a la plaque (W)
% Second membre
Fc = dt*mesh1.P1('x.^2+y.^2<0.2^2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%

%%% Avec laser
for i = 1:niter
    %%% Affichage
    mesh1.surf(U)
    caxis([0 0.04])
    pause(0.1)
    
    %%% Iteration
    U = Usolve(A, M, Fc, I, U);
end

% % On indique que le laser va etre coupe
% mesh1.surf(zeros(nP,1));
% pause(1);
% 
% %%% Sans laser
% for i = 11:niter
%     %%% Affichage
%     mesh1.surf(U)
%     caxis([0 0.04])
%     pause(0.1)
%     
%     %%% Iteration
%     U = Usolve(A, M, zeros(nP,1), I, U);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Affichage au temps final %%%%%%%%%
% mesh1.surf(U)
% caxis([0 1000])
% 
% figure;
% mesh1.plot

% figure;
% mesh1.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%