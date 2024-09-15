clear all 
a=0;
b=10;
n=1000;
h=(b-a)/n;
t =[a,a,a,a,a,a:h:b,b,b,b,b,b]; % exemple de n�uds
k =6; % exemple de degr�
x = linspace(t(k), t(end), 10000); % valeurs o� �valuer la B-spline
[B,K]= all_bsplines(x, k, t,h);
y=sin(x);
T=K*B;
plot(x,y,'r',x,T,'b')
E=y-T;
H=max(E)
L2_norm = norm(E,2)
L_infini = norm(E,inf)