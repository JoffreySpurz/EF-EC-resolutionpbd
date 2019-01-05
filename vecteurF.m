%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P : liste des sommets du maillage, pour la ieme ligne : Pi = [Pi1 Pi2]
% T : liste des triangles du maillage, pour la ieme ligne : Ti = [nPi1 nPi2
% nPi3]
% Flaser : vecteur laser sur les triangles
% dt : pas en espace
% F : vecteur second membre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ F ] = vecteurF( P, T, Flaser, dt )
% Nombre de triangles
NT = length(T);

% Nombre de sommets
NP = length(P);

% Innitialisation de F
F = zeros(NP, 1);

% Iteration sur les triangles
for i=1:NT
    % On recupere Pi
    P1 = P(T(i,1),:);
    P2 = P(T(i,2),:);
    P3 = P(T(i,3),:);
    
    % On definit dn et gn
    [S,gn] = geo(P1', P2', P3');
            
    % On recupere a, b
    [a, b] = base(P1', P2', P3');
    
    % On calcule Aik,il
    for k=1:3
            F(T(i,k)) = F(T(i,k)) + S*(a(1:2, k)' * gn(1:2) + b(1,k)) * Flaser(i) * dt;
    end
end

end

