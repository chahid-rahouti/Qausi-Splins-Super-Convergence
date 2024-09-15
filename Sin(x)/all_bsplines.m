function [B,K] = all_bsplines(x, k, t,h)
% x: vecteur de points d'évaluation des B-splines
% k: ordre des B-splines (k=4 pour les B-splines cubiques)
% t: vecteur de noeuds
for i = 1:length(t)-k
    K(1)=sin(t(6));
    %K(2)=(17/105)*sin((t(5)+t(6))/2)+(35/32)*sin((t(6)+t(7))/2)-(35/96)*sin((t(7)+t(8))/2)+(21/160)*sin((t(8)+t(9))/2)-(5/224)*sin((t(9)+t(10))/2);
    K(2)=(163/300)*sin(t(6))+sin(t(7))-sin(t(8))+(2/3)*sin(t(9))-(1/4)*sin(t(10))+(1/25)*sin(t(11));
    %K(3)=(-19/45)*sin((t(5)+t(6))/2)+(377/288)*sin((t(6)+t(7))/2)+(61/288)*sin((t(7)+t(8))/2)-(59/480)*sin((t(8)+t(9))/2)+(7/288)*sin((t(9)+t(10))/2);
    K(3)=(1/200)*sin(t(6))+(103/60)*sin(t(7))-(73/60)*sin(t(8))+(7/10)*sin(t(9))-(29/120)*sin(t(10))+(11/300)*sin(t(11));
    %K(4)=(47/315)*sin((t(5)+t(6))/2)-(77/144)*sin((t(6)+t(7))/2)+(251/144)*sin((t(7)+t(8))/2)-(97/240)*sin((t(8)+t(9))/2)+(47/1008)*sin((t(9)+t(10))/2);
    K(4)=(-41/400)*sin(t(6))+(43/60)*sin(t(7))+(103/120)*sin(t(8))-(7/10)*sin(t(9))+(13/48)*sin(t(10))-(13/300)*sin(t(11));
    K(i)= (13/240)*(sin(t(i+1))+sin(t(i+5)))-(7/15)*(sin(t(i+2))+sin(t(i+4)))+(73/40)*(sin(t(i+3)));
    K(length(t)-k-3)=(-41/400)*(sin(t(length(t)-k+1)))+(43/60)*(sin(t(length(t)-k)))+(103/120)*(sin((t(length(t)-k-1))))-(7/10)*(sin(t(length(t)-k-2)))+(13/48)*(sin(t(length(t)-k-3)))-(13/300)*(sin(t(length(t)-k-4)));
    K(length(t)-k-2)=(1/200)*(sin(t(length(t)-k+1)))+(103/60)*(sin(t(length(t)-k)))-(73/60)*(sin(t(length(t)-k-1)))+(7/10)*(sin((t(length(t)-k-2))))-(29/120)*(sin(t(length(t)-k-3)))+(11/300)*(sin(t(length(t)-k-4)));
    K(length(t)-k-1)=(163/300)*(sin(t(length(t)-k+1)))+(sin(t(length(t)-k)))-(sin(t(length(t)-k-1)))+(2/3)*sin((t(length(t)-k-2)))-(1/4)*sin(t(length(t)-k-3))+(1/25)*(sin(t(length(t)-k-4)));
    K(length(t)-k)=sin(t(length(t)-k+1));
end
% Initialisation de la matrice B
B = zeros(length(t)-k, length(x));

% Calcul de toutes les B-splines pour chaque point d'évaluation
for i = 1:length(t)-k
    B(i,:) = bspline(x, i, k, t,h);
end
end
function B = bspline(x, i, k, t,h)
% x: point d'évaluation de la B-spline
% i: indice de la B-spline
% k: ordre de la B-spline (k=4 pour les B-splines cubiques)
% t: vecteur de noeuds

% Initialisation du vecteur B
B = zeros(size(x));

% Cas de base
if k == 1
    % La B-spline est égale à 1 si x est compris entre les noeuds i et i+1
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
    %Calcul de la B-spline à partir des B-splines de l'ordre précédent
    alpha = (x - t(i)) ./ (t(i+k-1) - t(i));
    beta = (t(i+k) - x) ./ (t(i+k) - t(i+1));
    B = alpha .* bspline(x, i, k-1, t,h) + beta .* bspline(x, i+1, k-1, t,h);
    
    end
    end
end
end





