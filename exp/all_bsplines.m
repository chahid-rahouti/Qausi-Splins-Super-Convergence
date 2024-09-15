function [B,K] = all_bsplines(x, k, t,h)
% x: vecteur de points d'�valuation des B-splines
% k: ordre des B-splines (k=4 pour les B-splines cubiques)
% t: vecteur de noeuds
for i = 1:length(t)-k
    K(1) = (-1/2)*(exp(t(6)/2)-1);
    %K(2) = (17/105)*(-1/2)*(exp((t(5)+t(6))/2/2)-1) + (35/32)*(-1/2)*(exp((t(6)+t(7))/2/2)-1) - (35/96)*(-1/2)*(exp((t(7)+t(8))/2/2)-1) + (21/160)*(-1/2)*(exp((t(8)+t(9))/2/2)-1) - (5/224)*(-1/2)*(exp((t(9)+t(10))/2/2)-1);
    K(2) = (163/300)*(-1/2)*(exp(t(6)/2)-1) + (-1/2)*(exp(t(7)/2)-1) - (-1/2)*(exp(t(8)/2)-1) + (2/3)*(-1/2)*(exp(t(9)/2)-1) - (1/4)*(-1/2)*(exp(t(10)/2)-1) + (1/25)*(-1/2)*(exp(t(11)/2)-1);
    %K(3) = (-19/45)*(-1/2)*(exp((t(5)+t(6))/2/2)-1) + (377/288)*(-1/2)*(exp((t(6)+t(7))/2/2)-1) + (61/288)*(-1/2)*(exp((t(7)+t(8))/2/2)-1) - (59/480)*(-1/2)*(exp((t(8)+t(9))/2/2)-1) + (7/288)*(-1/2)*(exp((t(9)+t(10))/2/2)-1);
    K(3) = (1/200)*(-1/2)*(exp(t(6)/2)-1) + (103/60)*(-1/2)*(exp(t(7)/2)-1) - (73/60)*(-1/2)*(exp(t(8)/2)-1) + (7/10)*(-1/2)*(exp(t(9)/2)-1) - (29/120)*(-1/2)*(exp(t(10)/2)-1) + (11/300)*(-1/2)*(exp(t(11)/2)-1);
    %K(4) = (47/315)*(-1/2)*(exp((t(5)+t(6))/2/2)-1) - (77/144)*(-1/2)*(exp((t(6)+t(7))/2/2)-1) + (251/144)*(-1/2)*(exp((t(7)+t(8))/2/2)-1) - (97/240)*(-1/2)*(exp((t(8)+t(9))/2/2)-1) + (47/1008)*(-1/2)*(exp((t(9)+t(10))/2/2)-1);
    K(4) = (-41/400)*(-1/2)*(exp(t(6)/2)-1) + (43/60)*(-1/2)*(exp(t(7)/2)-1) + (103/120)*(-1/2)*(exp(t(8)/2)-1) - (7/10)*(-1/2)*(exp(t(9)/2)-1) + (13/48)*(-1/2)*(exp(t(10)/2)-1) - (13/300)*(-1/2)*(exp(t(11)/2)-1);
    K(i) = (13/240)*((-1/2)*(exp(t(i+1)/2)-1) + (-1/2)*(exp(t(i+5)/2)-1)) - (7/15)*((-1/2)*(exp(t(i+2)/2)-1) + (-1/2)*(exp(t(i+4)/2)-1)) + (73/40)*((-1/2)*(exp(t(i+3)/2)-1));
    K(length(t)-k-3) = (-41/400)*((-1/2)*(exp(t(length(t)-k+1)/2)-1)) + (43/60)*((-1/2)*(exp(t(length(t)-k)/2)-1)) + (103/120)*((-1/2)*(exp(t(length(t)-k-1)/2)-1)) - (7/10)*((-1/2)*(exp(t(length(t)-k-2)/2)-1)) + (13/48)*((-1/2)*(exp(t(length(t)-k-3)/2)-1)) - (13/300)*((-1/2)*(exp(t(length(t)-k-4)/2)-1));
    K(length(t)-k-2) = (1/200)*((-1/2)*(exp(t(length(t)-k+1)/2)-1)) + (103/60)*((-1/2)*(exp(t(length(t)-k)/2)-1)) - (73/60)*((-1/2)*(exp(t(length(t)-k-1)/2)-1)) + (7/10)*((-1/2)*(exp(t(length(t)-k-2)/2)-1)) - (29/120)*((-1/2)*(exp(t(length(t)-k-3)/2)-1)) + (11/300)*((-1/2)*(exp(t(length(t)-k-4)/2)-1));
    K(length(t)-k-1) = (163/300)*((-1/2)*(exp(t(length(t)-k+1)/2)-1)) + ((-1/2)*(exp(t(length(t)-k)/2)-1)) - ((-1/2)*(exp(t(length(t)-k-1)/2)-1)) + (2/3)*((-1/2)*(exp(t(length(t)-k-2)/2)-1)) - (1/4)*((-1/2)*(exp(t(length(t)-k-3)/2)-1)) + (1/25)*((-1/2)*(exp(t(length(t)-k-4)/2)-1));
    K(length(t)-k) = (-1/2)*(exp(t(length(t)-k+1)/2)-1);
end
% Initialisation de la matrice B
B = zeros(length(t)-k, length(x));

% Calcul de toutes les B-splines pour chaque point d'�valuation
for i = 1:length(t)-k
    B(i,:) = bspline(x, i, k, t,h);
end
end
function B = bspline(x, i, k, t,h)
% x: point d'�valuation de la B-spline
% i: indice de la B-spline
% k: ordre de la B-spline (k=4 pour les B-splines cubiques)
% t: vecteur de noeuds

% Initialisation du vecteur B
B = zeros(size(x));

% Cas de base
if k == 1
    % La B-spline est �gale � 1 si x est compris entre les noeuds i et i+1
    B(x >= t(i) & x <= t(i+1)) = 1; 
else if (t(i+k-1) == t(i))
    alpha = 0;
    beta = ((t(i+k) - x) ./ h) ;
    B = alpha .* bspline(x, i, k-1, t,h) + beta.* bspline(x, i+1, k-1, t,h);
else if ((t(i+k) == t(i+1)))
     alpha = ((x - t(i)) ./ h );
    beta = 0 ;
    B = alpha .* bspline(x, i, k-1, t,h) + beta .* bspline(x, i+1, k-1, t,h);
  else
    %Calcul de la B-spline � partir des B-splines de l'ordre pr�c�dent
    alpha = (x - t(i)) ./ (t(i+k-1) - t(i));
    beta = (t(i+k) - x) ./ (t(i+k) - t(i+1));
    B = alpha .* bspline(x, i, k-1, t,h) + beta .* bspline(x, i+1, k-1, t,h);
    
    end
    end
end
end





