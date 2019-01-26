clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Maillage %%%%%%%%%%
% Avec (1) ou sans (0) fissure ?
isfissure = 1;

% Fissure carree
pos_x = -0.5;
pos_y = 0.25;
longueur_y = 0.1;
longueur_x = 1;

% 1 Unit square
if isfissure == 1
    nodes = [-1 -1; 1 -1 ; 1 1; -1 1;pos_x pos_y ; pos_x+longueur_x pos_y; pos_x+longueur_x pos_y+longueur_y ;pos_x pos_y+longueur_y];
    edges = {{1,2} {2,3} {3,4} {4,1} {5,6,'lr'} {6,7,'lr'} {7,8,'lr'} {8,5,'lr'}}; 
else
    nodes = [-1 -1; 1 -1 ; 1 1; -1 1];
    edges = {{1,2} {2,3} {3,4} {4,1}};
end
domain1 = Domain(nodes,edges);

% mesh
dx = 0.05; % taille d'un element
mesh1 = Mesh(domain1,dx,'save');

% Liste des sommets
P = mesh1.nodes;

% Liste des triangles
T = mesh1.triangles;

% Liste des points du bords
I = mesh1.boundary([1 2 3 4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Parametres %%%%%%%%%%
% Nombre de sommets
nP = length(P); 

% Nombre de triangles
nT = length(T);

% Changements d'echelles
l = 10^(-2);% Convertion des m en cm

% Conductivite thermique (W/m/K)
lambda_plaque = 50.2; % valeur de l'acier 
lambda_fiss = 0.024; %  valeur de l'air

% Masse volumique (Kg/m^3)
rho_plaque = 7850;% valeur de l'acier
rho_fiss = 1;% valeur de l'air

% Capacite thermique massique (J/K/Kg)
cp_plaque = 444;% valeur de l'acier
cp_fiss = 1004;% valeur de l'air

% Beta = Rho*cp (J/K/m^3)
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
beta_fiss = rho_fiss*cp_fiss ; % rho*cp dans la fissure

% Plaque avec fissures
if isfissure == 1
    alpha = 1* mesh1.P0([ num2str(pos_x) '<=x']).*mesh1.P0(['x<=' num2str(pos_x+longueur_x)]).*mesh1.P0([num2str(pos_y) '<=y']).*mesh1.P0(['y<= ' num2str(pos_y+longueur_y)]); % Indicatrice : 1 dans la fissure 0 sinon
else
    alpha = mesh1.P0(0);
end

Beta = (beta_plaque/lambda_plaque)*(1+(beta_fiss-beta_plaque)/beta_plaque*alpha); % Vecteur Rho*cp sur chaque triangle
Lambda = 1+(lambda_fiss/lambda_plaque-1)*alpha; % Vecteur Conductivite thermique (lambda) sur les triangles

% Iteration
T = 10;% Temps (s) final d'experience
Tstop = floor(T/3);% Temps (s) d'arret du laser
dt =T/100; % Pas en temps
niter = floor(T/dt); % Nombre d'iterations pour la resolution
niter_laser = floor(Tstop/dt); % Nombre d'iterations pour le laser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Matrices d'iterations %%%%%%%%
%%% Raideur
% Paramètre : Rho*cp sur chaque triangle
LAMBDA = mesh1.P0(Lambda,2);
% Matrice de raideur
Kc = dt*mesh1.stiffness(LAMBDA);%matriceK(P, T, Lambda, dt);


%%% Masse
% Paramètre : Conductivite thermique (lambda) sur les triangles
BETA = mesh1.P0(Beta*l^2);
% Matrice de masse
Mc = mesh1.mass(BETA);
M = mesh1.mass();
%mesh1.mass(Beta.*ones(nT,1));


%%% Membre de gauche
% Probleme AX(n) = McX(n-1) + Fc
A = Mc + Kc;

%%% Prise en compte des CL %%%
% Nombre de points du bord
nI = length(I);

% CL de Dirichlet
A(I,:) = 0;
A(I,I) = speye(nI);
M(I,:) = 0;
M(I,I) = speye(nI);
Mc(I,:) = 0;
Mc(I,I) = speye(nI);

%disp('Construction de B')

% Probleme X(n) = A\M ( X(n-1) + F ) = B ( X(n-1) + F )
%B = A\M;


%%% Laser
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv*l^3/(lambda_plaque*l); % Puissance volumique du laser utilisee
% Second membre 
sigma = 0.005;
Fc = M*dt*S/sqrt(2*pi*sigma^2)*mesh1.P1(['exp(-0.5*(x.^2+y.^2)/' num2str(sigma^2) ')']);% Gaussien
%Fc=dt*S*mesh1.P1('x.^2+y.^2<0.2^2');% Cercle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Debut des iterations') 

%%%%%%%%% Resolution %%%%%%%%%%%
%%% Avec laser
U = Usolve(A, Mc, Fc, nP, niter, niter_laser);

%%% Save
if isfissure ==1
    save('Probleme_direct_avec_fissure_rectangulaire','U', 'mesh1', 'niter_laser')
else 
    save('Probleme_direct_sans_fissure','U', 'mesh1', 'niter_laser')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%