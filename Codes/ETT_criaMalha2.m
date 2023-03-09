function [ Malha ] = ETT_criaMalha2( geom , ps )
%% Carregando Valores %%
% ================================ AVIÃO ================================ %
b = geom.b;
b_sec = geom.b_sec;
perfil = geom.perfil;
i_w = geom.i_w;
xBA = geom.xBA;
zBA = geom.zBA;

spline_c = geom.spline_c;
spline_l = geom.spline_l;
spline_tw = geom.spline_tw;

%% Corpo da Função %%
% ======================= Criando geometria base ======================== %
n = 100;
m = 100;
y = linspace( 0 , b/2 , n );
x = linspace( 0 , 1 , m )';

Zs = zeros(m,n);
Zi = zeros(m,n);
Zs_sup = zeros(m,n);
Zi_sup = zeros(m,n);

fb = b_sec / b;
num_y = round( n * fb );
num_y(end) = n - sum( num_y(1:end-1) );
num_y = [ 1 cumsum(num_y) ];

for i = 1:length(b_sec)
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );
    
    num_pontos_sec = num_y(i+1) - num_y(i) + 1;
    
    y( num_y(i):num_y(i+1) ) = linspace( y0 , y1 , num_pontos_sec );
    
    Z0s = c0 * interp1( ps(perfil(i)).pontos_sup(:,1) , ps(perfil(i)).pontos_sup(:,2) , x , 'linear' , 'extrap' );
    Z0i = c0 * interp1( ps(perfil(i)).pontos_inf(:,1) , ps(perfil(i)).pontos_inf(:,2) , x , 'linear' , 'extrap' );
    Z1s = c1 * interp1( ps(perfil(i+1)).pontos_sup(:,1) , ps(perfil(i+1)).pontos_sup(:,2) , x , 'linear' , 'extrap' );
    Z1i = c1 * interp1( ps(perfil(i+1)).pontos_inf(:,1) , ps(perfil(i+1)).pontos_inf(:,2) , x , 'linear' , 'extrap' );
    
    
    Zs(:,num_y(i):num_y(i+1)) = linspaceVet( Z0s , Z1s , num_pontos_sec );
    Zi(:,num_y(i):num_y(i+1)) = linspaceVet( Z0i , Z1i , num_pontos_sec );
    
    %
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );
    
    
end

x_14 = ppval(spline_c,0)/4 + [ 0 , cumsum(tand( ppval(spline_l,y(1:end-1)) ) .* diff( y ))]';

X = linspaceVet( x_14 - ppval(spline_c,y')/4 , x_14 + 3*ppval(spline_c,y')/4 , m )';
Y = repmat( y , m , 1 );

% =============================== Twists ================================ %
Z_med = (Zs + Zi)/2;

for i = 1:size(Z_med,2)
   
   centro_twist = [interp1( linspace(0,1,m) , X(:,i) , 0 ) , interp1( linspace(0,1,m) , Z_med(:,i) , 0 )];
   [ x_twist , z_sup_twist , z_inf_twist ] = aplicaTwist( X(:,i)' , Zs(:,i)' , Zi(:,i)' , i_w , centro_twist );
   
   X(:,i) = x_twist;
   Zi(:,i) = z_inf_twist;
   Zs(:,i) = z_sup_twist;
   
   twist = ppval( spline_tw , Y(1,i) );
   
   centro_twist = [interp1( linspace(0,1,m) , X(:,i) , 1/4 ) , interp1( linspace(0,1,m) , Z_med(:,i) , 1/4 )];
   [ x_twist , z_sup_twist , z_inf_twist ] = aplicaTwist( X(:,i)' , Zs(:,i)' , Zi(:,i)' , twist , centro_twist );
   
   X(:,i) = x_twist;
   Zi(:,i) = z_inf_twist;
   Zs(:,i) = z_sup_twist;
        
    
end

% =============================== Alturas =============================== %
Zi = Zi + zBA;
Zs = Zs + zBA;
Zm = (Zi + Zs)/2;

% ================================== X ================================== %
X = X + xBA;

%% ATUALIZANDO STRUCT %%
Malha(1).X_orig = X;
Malha(1).Z_sup_orig = Zs;
Malha(1).Z_inf_orig = Zi;
Malha(1).Z_med_orig = Zm;

Malha(1).X = X;
Malha(1).Y = Y;
Malha(1).Z_sup = Zs;
Malha(1).Z_inf = Zi;
Malha(1).Z_med = Zm;

% Malha(1).posMalha = posMalha;


end


%% Função linspaceVet() %%
function y = linspaceVet(d1,d2,n)
n1 = n-1;
N = repmat(0:n1,length(d1),1);
y = d1 + N.*(d2-d1)./n1;
end

%% Função rgb() %%
function [rgb_norm] = rgb( r , g , b )

rgb_norm = [r , g , b] /256;

end

%% Função aplicaTwist() %%
% ============================== DESCRIÇÃO ============================== %
% Essa função recebe dois vetores referentes às coordenadas X e Z do perfil
% e aplica um twist.
function [ x_twist , z_sup_twist , z_inf_twist ] = aplicaTwist( x , z_sup , z_inf , twist , centro_twist )
x = x - centro_twist(1);
z_sup = z_sup - centro_twist(2);
z_inf = z_inf - centro_twist(2);

R = [cosd(-twist) -sind(-twist); sind(-twist) cosd(-twist)];
v_sup = [x;z_sup];
v_inf = [x;z_inf];

v_sup_twist = R*v_sup;
v_inf_twist = R*v_inf;

x_twist = v_sup_twist(1,:);
z_sup_twist = v_sup_twist(2,:);
z_inf_twist = interp1(v_inf_twist(1,:),v_inf_twist(2,:),x_twist,'linear','extrap');

x_twist = x_twist + centro_twist(1);
z_sup_twist = z_sup_twist + centro_twist(2);
z_inf_twist = z_inf_twist + centro_twist(2);

end