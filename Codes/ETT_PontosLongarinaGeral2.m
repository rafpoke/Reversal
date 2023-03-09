function [ L ] = ETT_PontosLongarinaGeral2( geom , L, ps )
%% CARREGANDO VALORES %%
Malha = {geom.Malha};
is_long = [L.is_long];
% is_trelica = [L.is_trelica];
config = {L(is_long).config};
% fcC = {Longarina.fcC};
fc = {L(is_long).fc};
fb = {L(is_long).fb};
fb0 = [L(is_long).fb0];
% fhC = {Longarina.fhC};
fh = {L(is_long).fh};
S_fc = {L(is_long).S_fc};
num_pontos = [L(is_long).num_pontos];
idx_malha = [L(is_long).idx_malha];
l_long = {L(is_long).l_long};
d_long = {L(is_long).d_long};


num_long = length(config);
num_tot = length(L);

% [ L ] = ETT_criaTrelicado( L(is_trelica) );

%% CORPO DA FUNÇÃO %%
aux = 1:num_long-5;
aux2 = 1:num_tot;
aux2 = aux2(is_long);
% aux = aux(is_long);


% Para simplificar os programas, cada longarina é tratada separadamente
for i = aux
    k = aux2(i);
    L(k).pontos = zeros( num_pontos(i) , 3 );
    num_sec = length( L(k).config );

    [ L(k).pontos_sec , L(k).posLong , L(k).Hmax_sec ] = pontosLongarinaSec( Malha{idx_malha(i)} , fb{i} , fc{i} , fh{i} , num_pontos(i) , fb0(i) , l_long{i} , d_long{i} );

    % calculo de enflechamentos
    l = nan(1,num_sec+1);
    l(1) = 0;
    l(end) = calcula_l_end( Malha{idx_malha(i)} , fc{i} );
    is_perfil = strcmp( config{i}, 'BA' );
    if is_perfil
            fBA = L(k).fBA;
            nx = L(k).nx;
            L(k).Malha_Long = ETT_criaMalhaBA2( geom(idx_malha(i)) , fBA, fb{i}, num_pontos(i), nx , ps );
            L(k).pontos = Acha_Centroide(L(k).Malha_Long);
    else
        for j = 1:num_sec
            num_pontos_sec = diff( L(k).posLong(j:j+1) ) + 1;
            idx_pontos = L(k).posLong(j):L(k).posLong(j+1);
            [ L(k).pontos(idx_pontos,:) , l(j+1) ] = pontosLongarinaSimples( config{i}{j} , L(k).pontos_sec(j:j+1,:) , num_pontos_sec , S_fc{i}(j) , l(j) , l(j+1) );
        end
    end
end
end

%% ETT_criaTrelicado %%
function [L] = ETT_criaTrelicado(L)
P0 = L.P0;
P4 = L.P4;
f1 = L.f1;
f2 = L.f2;
f3 = L.f3;
num_T = L.num_T;
b = L.b;
h = L.h;
Sh = L.Sh;
Sb = L.Sb;
num_pontos = L.num_pontos;

config = L.config;

%%  Diametros possíveis %%
D_table = ...
[2 1
 3 2
 4 3
 4 2
 5 3
 6 4
 8 6]*1e-3;

%% Cria sp %%
sp = struct;

sp.config = config;

sp.b = b;
sp.Sb = Sb;
sp.h = h;
sp.Sh = Sh;
sp.P0 = P0;
sp.P4 = P4;
sp.f1 = f1;
sp.f2 = f2;
sp.f3 = f3;
sp.num_T = num_T;
sp.num_pontos = num_pontos;

sp.material = {materiais.carbono};

%% Configuração %%
num = 28 * sp.num_T + 2;

sp.cc = cell(1,num); 
sp.connect = cell(1,num);

if strcmp(config,'plano')
    idx_D = ceil( rand(1,5*num_T) * 7 );
    idx_D(:) = 7;
    bool = true(1,5*num_T);
    for i = 1:4
        sp.cc{1 + num_T*(i-1)} = { 1 + num_T*(i-1) , 1 , 1:6 };
    end

elseif strcmp(config,'triangulo')
    idx_D = ceil( rand(1,12*num_T) * 5);
    idx_D(:) = 4;
    bool = true(1,12*num_T);
    for i = 1:9
        sp.cc{1 + num_T*(i-1)} = { 1 + num_T*(i-1) , 1 , 1:6 };
    end
    
elseif strcmp(config,'quadrado')
    idx_D = ceil( rand(1,22*num_T) * 7 );
    idx_D(:) = 3;
    bool = true(1,22*num_T);
    for i = 1:22
        sp.cc{1 + num_T*(i-1)} = { 1 + num_T*(i-1) , 1 , 1:6 };
%         sp.cc{num_T + num_T*(i-1)} = { num_T + num_T*(i-1) , 1 , 1:6 };
    end
end

D_ext = D_table(idx_D,1);
D_int = D_table(idx_D,2);

sp.D_ext = D_ext;
sp.D_int = D_int;


sp.bool = logical( bool );

L.sp = sp;

end

%% Acha_Centriode %%
function [pontos] = Acha_Centroide(Malha_Long)
[m,n] = size(Malha_Long.Y);

xc = zeros(n,1);
yc = (Malha_Long.Y(1,:))';
zc = zeros(n,1);

for i = 1:n
    A = 0;
    CX = 0;
    CZ = 0;
    X = [flip(Malha_Long.X(:,i));Malha_Long.X(2:end,i);Malha_Long.X(end,i)]';
    Z = [flip(Malha_Long.Z_sup(:,i));Malha_Long.Z_inf(2:end,i);Malha_Long.Z_sup(end,i)]';
    for j = 1:2*m-1
        A = 0.5*(X(j)*Z(j+1) - X(j+1)*Z(j)) + A;
    end
    for j = 1:2*m-1
        CX = (X(j) + X(j+1))*(X(j)*Z(j+1) - X(j+1)*Z(j)) + CX;
        CZ = (Z(j) + Z(j+1))*(X(j)*Z(j+1) - X(j+1)*Z(j)) + CZ;
    end
    xc(i) = CX/(6*A);
    zc(i) = CZ/(6*A);
end
pontos = [xc yc zc];
end
%% ETT_criaMalhaBA %%
function [ Malha ] = ETT_criaMalhaBA2( geom_Asa , fBA, fb, n, m , ps )
%% Carregando Valores %%
b = geom_Asa.b;
b_sec = geom_Asa.b_sec;
perfil = geom_Asa.perfil;
i_w = geom_Asa.i_w;
xBA = geom_Asa.xBA;
zBA = geom_Asa.zBA;
spline_c = geom_Asa.spline_c;
spline_l = geom_Asa.spline_l;
spline_tw = geom_Asa.spline_tw;
% 
% fBA = geom_BA.fBA;
% fb = geom_BA.fb;
% n = geom_BA.n;
% m = geom_BA.m;


%% Corpo da Função %%
% ======================= Criando geometria base ======================== %
n_sec_long = length(fBA)-1;

b_sec_long = b*fb;
b_sec_abs = cumsum(b_sec);

b_sec_novo = sort([b_sec_abs,b_sec_long(1:end-1)]);

sec_long = zeros(1,n_sec_long-1);
for i = 1:n_sec_long-1
    sec_long(i) = find(b_sec_novo == b_sec_long(i));
    perfil = [perfil(1:sec_long(i)), perfil(sec_long(i)), perfil(sec_long(i)+1:end)];
end
% 
% b_sec_novo = diff(b_sec_novo);
% b_sec_novo = [b_sec_novo,b-sum(b_sec_novo)];

sec_long = [sec_long,length(b_sec_novo)];

y = linspace( 0 , b/2 , n );
% vet_fBA = linspace(fBA(1),fBA(2),n);

% x = linspaceVet(zeros(n_sec,1),fBA',m)';
% x = linspace( 0 , 1 , m )';

% Zs = zeros(m,n);
% Zi = zeros(m,n);
Zs = [];
Zi = [];
Zs_sup = zeros(m,n);
Zi_sup = zeros(m,n);

fb = b_sec_novo / b;
num_y = round( n * fb );
interv = num_y;

% num_y(end) = n - sum( num_y(1:end-1) );
% num_y = [ 1 cumsum(num_y) ];

num_y = [1 num_y];

for i = 1:length(b_sec_novo)
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );
    
    num_pontos_sec = num_y(i+1) - num_y(i) + 1;
    
    y( num_y(i):num_y(i+1) ) = linspace( y0 , y1 , num_pontos_sec );
end

x_14 = ppval(spline_c,0)/4 + [ 0 , cumsum(tand( ppval(spline_l,y(1:end-1)) ) .* diff( y ))]';

corda = ppval(spline_c,y');

BFBA = [];
raiz_long = fBA(1)*corda(1);
inter_frac = 0;
for i = 1:n_sec_long
    BFBA = [BFBA, linspace(raiz_long,fBA(i+1)*corda(num_y(sec_long(i))),abs(inter_frac - interv(sec_long(i))))];
    raiz_long = fBA(i+1)*corda(num_y(sec_long(i)));
    inter_frac = interv(sec_long(i));
end
BFBA = BFBA';
x_BA = x_14 - corda/4;
x_BF = x_BA + BFBA;
X = linspaceVet( x_BA , x_BF , m )';

frac_x = (X(end,:)' - X(1,:)')./corda;
frac_x = frac_x(num_y);

Y = repmat( y , m , 1 );

for i = 1:length(b_sec_novo)    
    
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );  
    
    num_pontos_sec = num_y(i+1) - num_y(i) + 1;
    
    x0 = linspace(0,frac_x(i),m)';
    x1 = linspace(0,frac_x(i+1),m)';
    Z0s = c0 * interp1( ps(perfil(i)).pontos_sup(:,1) , ps(perfil(i)).pontos_sup(:,2) , x0 , 'linear' , 'extrap' );
    Z0i = c0 * interp1( ps(perfil(i)).pontos_inf(:,1) , ps(perfil(i)).pontos_inf(:,2) , x0 , 'linear' , 'extrap' );
    Z1s = c1 * interp1( ps(perfil(i+1)).pontos_sup(:,1) , ps(perfil(i+1)).pontos_sup(:,2) , x1 , 'linear' , 'extrap' );
    Z1i = c1 * interp1( ps(perfil(i+1)).pontos_inf(:,1) , ps(perfil(i+1)).pontos_inf(:,2) , x1 , 'linear' , 'extrap' );
    
    
%     Zs = [Zs;linspaceVet( Z0s , Z1s , num_pontos_sec )];
%     Zi = [Zi;linspaceVet( Z0i , Z1i , num_pontos_sec )];
    Zs(:,num_y(i):num_y(i+1)) = linspaceVet( Z0s , Z1s , num_pontos_sec );
    Zi(:,num_y(i):num_y(i+1)) = linspaceVet( Z0i , Z1i , num_pontos_sec );
    
    %
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );
    
    
end

%% Asa %%
Zs_Asa = zeros(m,n);
Zi_Asa = zeros(m,n);
x = linspace( 0 , 1 , m )';

for i = 1:length(b_sec_novo)
    y0 = (num_y(i)/n)*b/2 * (i~=1);
    y1 = (num_y(i+1)/n)*b/2;
    c0 = ppval( spline_c , y0 );
    c1 = ppval( spline_c , y1 );
    
    num_pontos_sec = num_y(i+1) - num_y(i) + 1;
    
    y( num_y(i):num_y(i+1) ) = linspace( y0 , y1 , num_pontos_sec );
    
    Z0s_Asa = c0 * interp1( ps(perfil(i)).pontos_sup(:,1) , ps(perfil(i)).pontos_sup(:,2) , x , 'linear' , 'extrap' );
    Z0i_Asa = c0 * interp1( ps(perfil(i)).pontos_inf(:,1) , ps(perfil(i)).pontos_inf(:,2) , x , 'linear' , 'extrap' );
    Z1s_Asa = c1 * interp1( ps(perfil(i+1)).pontos_sup(:,1) , ps(perfil(i+1)).pontos_sup(:,2) , x , 'linear' , 'extrap' );
    Z1i_Asa = c1 * interp1( ps(perfil(i+1)).pontos_inf(:,1) , ps(perfil(i+1)).pontos_inf(:,2) , x , 'linear' , 'extrap' );
    
    
    Zs_Asa(:,num_y(i):num_y(i+1)) = linspaceVet( Z0s_Asa , Z1s_Asa , num_pontos_sec );
    Zi_Asa(:,num_y(i):num_y(i+1)) = linspaceVet( Z0i_Asa , Z1i_Asa , num_pontos_sec );
    
    
end

% =============================== Twists ================================ %
Z_med = (Zs + Zi)/2;
Z_med_Asa = (Zs_Asa + Zi_Asa)/2;
X_Asa = linspaceVet( x_14 - ppval(spline_c,y')/4 , x_14 + 3*ppval(spline_c,y')/4 , m )';


for i = 1:size(Z_med,2)
   
   centro_twist = [interp1( linspace(0,1,m) , X_Asa(:,i) , 0 ) , interp1( linspace(0,1,m) , Z_med_Asa(:,i) , 0 )];
   [ x_twist , z_sup_twist , z_inf_twist ] = aplicaTwist( X(:,i)' , Zs(:,i)' , Zi(:,i)' , i_w , centro_twist );
   
   X(:,i) = x_twist;
   Zi(:,i) = z_inf_twist;
   Zs(:,i) = z_sup_twist;
   
   twist = ppval( spline_tw , Y(1,i) );
   
   centro_twist = [interp1( linspace(0,1,m) , X_Asa(:,i) , 1/4 ) , interp1( linspace(0,1,m) , Z_med_Asa(:,i) , 1/4 )];
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
Malha.X_orig = X;
Malha.Z_sup_orig = Zs;
Malha.Z_inf_orig = Zi;
Malha.Z_med_orig = Zm;

Malha.X = X;
Malha.Y = Y;
Malha.Z_sup = Zs;
Malha.Z_inf = Zi;
Malha.Z_med = Zm;

% Malha(1).posMalha = posMalha;


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
%% Função pontosLongarinaSimples() %%

function [ pontos , l_end ] = pontosLongarinaSimples( config , pontos_sec , num_pontos , S_fc , l0 , l1 )
is_continuo = strcmp( config , 'sanduiche' );

if ~is_continuo
    % Longarina Discreta
    pontos = linspaceVet( pontos_sec(1,:)' , pontos_sec(2,:)' , num_pontos )';
    l_end = 0;

else
    % Longarina Contínua
    x0 = pontos_sec(1,1);
    x1 = pontos_sec(2,1);
    b_sec = diff( pontos_sec(:,2) );                                        % envergadura da seção de longarina
    
    if isnan(l1)
        pc = F_Spline_posx_mid( x0 , x1 ,  l0 , S_fc , b_sec );
    else
        pc = F_Spline_posx_end( x0 , x1 , l0 , l1 , b_sec );
    end
    
    x_long = polyval(pc,linspace(0,1,num_pontos));
    
    dx = x_long(end) - x_long(end-1);
    dy = b_sec / (num_pontos - 1);
    l_end = dx/dy;
    
    ph = F_Spline_posx(pontos_sec(1,3),pontos_sec(2,3),0,0);
    z_long = polyval(ph,linspace(0,1,num_pontos));
    
    y_long = linspace( pontos_sec(1,2) , pontos_sec(2,2) , num_pontos );
    
    pontos = [ x_long' , y_long' , z_long' ];
    
end

end

%% Função pontosLongarinaSec() %%
% ============================== DESCRIÇÃO ============================== %
% Essa função recebe uma malha de pontos e parâmetros geométricos da
% longarina e retorna uma matriz de pontos da longarina.

function [ pontos_sec , posLong , Hmax_sec] = pontosLongarinaSec( Malha , fb , fc , fh , num_pontos , fb0 , l_long , d_long )
%% CARREGANDO VALORES %%
X = Malha.X;
Y = Malha.Y;
Z_sup = Malha.Z_sup;
Z_inf = Malha.Z_inf;
Z_med = Malha.Z_med;

b = Y(1,end);                                                               % semi-envergadura da malha
espessura_perfil = Z_sup - Z_inf;                                           % espessura da malha
num_secoes = length(fb);


%% CORPO DA FUNÇÃO %%
% ======================== Pontos com Y Negativo ======================== %
if fb0 < 0 || fb(1) < 0
    flag_neg = 1;
    fb0 = -fb0;
    fb = -fb;
else
    flag_neg = 0;
end

% ============================ INICIALIZAÇÃO ============================ %
y_long = zeros(1,num_secoes+1);
y_long(1) = b * fb0;

% ======================== Posição em Y discreta ======================== %
for i = 2:num_secoes+1
    y_long(i) = min([ y_long(i-1) + (b - y_long(i-1)) * fb(i-1) , Y(1,end) ]);
end

% ====================== Posição em Y da longarina ====================== %
y_long_cont = linspace(y_long(1),y_long(end),num_pontos);

% =============================== posLong =============================== %
[~,posLong] = min( abs( y_long_cont' - y_long ) );                          % posLong é a posição em y mais próxima de y_disc (pode dar pau com seçoes muito pequenas)

% =================== CORDA NAS POSIÇÕES DA LONGARINA =================== %
cordas = interp1(Y(1,:),X(end,:)-X(1,:),y_long);

% ====================== POSIÇÃO EM X DA LONGARINA ====================== %
x_long = interp1(Y(1,:),X(1,:),y_long) + cordas .* fc;
idx = find(~isnan(l_long));
x_long(idx+1) = x_long(idx) + (y_long(idx+1)-y_long(idx)) .* tand(l_long(idx));

% ====================== POSIÇÃO EM Z DA LONGARINA ====================== %
z_long = griddata(X,Y,Z_med,x_long,y_long) + griddata(X,Y,espessura_perfil,x_long,y_long) .* fh;
idx = find(~isnan(d_long));
z_long(idx+1) = z_long(idx) + (y_long(idx+1)-y_long(idx)) .* tand(d_long(idx));

Hmax_sec = griddata(X,Y,espessura_perfil,x_long,y_long);
% ========================== Pontos Discretos =========================== %
pontos_sec = [ x_long' , y_long', z_long' ];

% ======================== Pontos com Y Negativo ======================== %
if flag_neg
    pontos_sec(:,2) = -pontos_sec(:,2);
end

end

%% Função calcula_l_end() %%
function [ l_end ] = calcula_l_end( Malha , fc )

l_end = (Malha.X(:,end) - Malha.X(:,end-1)) ./ (Malha.Y(:,end) - Malha.Y(:,end-1));
l_end = interp1(linspace(0,1,length(l_end)),l_end,fc(end));

end

%% Função linspaceVet() %%
function y = linspaceVet(d1,d2,n)
n1 = n-1;
N = repmat(0:n1,length(d1),1);
y = d1 + N.*(d2-d1)./n1;
end