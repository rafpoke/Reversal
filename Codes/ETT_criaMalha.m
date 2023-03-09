function [ Malha ] = ETT_criaMalha( geom , ps )
%% CARREGANDO VALORES %%
% ================================ AVIÃO ================================ %
b_total = geom.b;                                                           %[m] envergadura total da asa
b = geom.m_b;                                                               %[m] vetor de semi-envergaduras
c = geom.c_sec;                                                             %[m] vetor de cordas
l_BA = geom.l_BA;                                                           %[º] vetor de enflechamento no BA
perfil = geom.perfil;                                                       %[-] vetor de perfis
i_w1 = geom.i_w;                                                            %[º] incidência da aeronave
tw = geom.tw;
h_asa1_BA = geom.h_BA;
x_BA = geom.x_BA;

% ======================== CÁLCULOS GEOMÉTRICOS ========================= %
offset_BA_x = offsetBA_x(b,l_BA);                                           %[m] vetor offset do BA em x
cum_b = cumsum([0 b]);                                                      %[m] vetor com a posição em y das transições de seção
fb = 2*b / b_total;                                                         %[-] vetor com a fração de envergadura de cada seção

%% MALHA %%
twist = [0 cumsum(tw)];
% ========================= NÚMERO DE ELEMENTOS ========================= %
num_paineis_X = 20;                                                         % número de painéis em X
num_paineis_Y = 15;                                                         % número de painéis em Y
paineis_Y = round(fb*num_paineis_Y);                                        % painéis em Y em cada seção
cum_paineis_Y = [1 cumsum(paineis_Y)];
num_paineis_Y = sum(paineis_Y);
num_paineis_total = num_paineis_Y * num_paineis_X;
% ============================ INICIALIZAÇÃO ============================ %
posMalha = [1];
X = zeros(num_paineis_X,num_paineis_Y);
Z_sup = X;
Z_inf = X;
% ============================ DISTRIBUIÇÃO X =========================== %
theta = linspace(0,pi,num_paineis_X);
x = 0.5 * (1 - cos(theta));                                                 % distribuição cossenoidal ao longo da corda
% ================================ MALHA ================================ %
Y = repmat(linspace(0,cum_b(end),num_paineis_Y),num_paineis_X,1);
for secao = 1:length(paineis_Y)
    num_paineis = length(cum_paineis_Y(secao):cum_paineis_Y(secao+1));
%     X(:,cum_paineis_Y(secao):cum_paineis_Y(secao+1)) = linspaceVet(x'*c(secao)+offset_BA_x(secao),x'*c(secao+1)+offset_BA_x(secao+1),num_paineis);
    x0 = x*c(secao)+offset_BA_x(secao);
    x1 = x*c(secao+1)+offset_BA_x(secao+1);
    %ISSO AQUI INVERTE O PERFIL , MAIS FACIL QUE A FUNCAO INVERTEPERFIL QUE
    %EU PERDI TEMPO ESCREVENDO E ELA NEM FUNCIONA
    z0_sup = -c(secao)*interp1(ps(perfil(secao)).pontos_sup(:,1),ps(perfil(secao)).pontos_sup(:,2),x,'linear','extrap');
    z0_inf = -c(secao)*interp1(ps(perfil(secao)).pontos_inf(:,1),ps(perfil(secao)).pontos_inf(:,2),x,'linear','extrap');
    z1_sup = -c(secao+1)*interp1(ps(perfil(secao+1)).pontos_sup(:,1),ps(perfil(secao+1)).pontos_sup(:,2),x,'linear','extrap');
    z1_inf = -c(secao+1)*interp1(ps(perfil(secao+1)).pontos_inf(:,1),ps(perfil(secao+1)).pontos_inf(:,2),x,'linear','extrap');
    % ========================== APLICA TWIST =========================== %
    z_med = interp1(ps(perfil(secao)).pontos_sup(:,1),ps(perfil(secao)).pontos_sup(:,2),x,'linear','extrap')/2 + interp1(ps(perfil(secao)).pontos_inf(:,1),ps(perfil(secao)).pontos_inf(:,2),x,'linear','extrap')/2;
    z_med_quarto = interp1(x,z_med,1/4,'linear','extrap');
    
    centro_twist0 = [x0(1) + c(secao)/4 , c(secao)*z_med_quarto];
    centro_twist1 = [x1(1) + c(secao+1)/4 , c(secao+1)*z_med_quarto];
    [x0,z0_sup,z0_inf] = aplicaTwist( x0 , z0_sup , z0_inf , twist(secao) , centro_twist0 );
    [x1,z1_sup,z1_inf] = aplicaTwist( x1 , z1_sup , z1_inf , twist(secao+1) , centro_twist1 );
    X(:,cum_paineis_Y(secao):cum_paineis_Y(secao+1)) = linspaceVet(x0',x1',num_paineis);

    Z_sup(:,cum_paineis_Y(secao):cum_paineis_Y(secao+1)) = linspaceVet(z0_sup',z1_sup',num_paineis);
    Z_inf(:,cum_paineis_Y(secao):cum_paineis_Y(secao+1)) = linspaceVet(z0_inf',z1_inf',num_paineis);
    posMalha = [posMalha posMalha(end)+num_paineis-1];
end
% % ============================== SIMETRIA =============================== %
%  X = [fliplr(X) X(:,2:end)];
%  Y = [-fliplr(Y) Y(:,2:end)];
%  Z_sup = [fliplr(Z_sup) Z_sup(:,2:end)];
%  Z_inf = [fliplr(Z_inf) Z_inf(:,2:end)];

% ============================= INCIDÊNCIA ============================== %
for i = 1:length(Y(1,:))
    [x,z_sup,z_inf] = aplicaTwist(X(:,i)',Z_sup(:,i)',Z_inf(:,i)',i_w1,[0,0]);
    X(:,i) = x';
    Z_sup(:,i) = z_sup';
    Z_inf(:,i) = z_inf';
end
% =========================== SALVANDO DADOS ============================ %
X1_orig = X;
Y1 = Y;
%Perfil  invertido
 Z_sup1_orig = Z_inf;
 Z_inf1_orig = Z_sup;

%Perfil não invertido
%Z_sup1_orig = Z_sup;
%Z_inf1_orig = Z_inf;

% =========================== OFFSET DA MALHA =========================== %
espessura = 5e-3;                                                          % distância mínima entre a longarina e a asa
[ X1 , Z_sup1 , Z_inf1 ] = offsetMalha( X1_orig , Z_sup1_orig , Z_inf1_orig , espessura );

% % ============================= TRANSLAÇÕES ============================= %
Z_sup1 = Z_sup1 + h_asa1_BA;
Z_inf1 = Z_inf1 + h_asa1_BA;
Z_sup1_orig = Z_sup1_orig + h_asa1_BA;
Z_inf1_orig = Z_inf1_orig + h_asa1_BA;

X1 = X1 + x_BA;
X1_orig = X1_orig + x_BA;

%% LINHA MÉDIA %%
Z1_med_orig = (Z_sup1_orig + Z_inf1_orig)/2;

Z1_med = (Z_sup1 + Z_inf1)/2;

%% ATUALIZANDO STRUCT %%
Malha(1).X_orig = X1_orig;
Malha(1).Z_sup_orig = Z_sup1_orig;
Malha(1).Z_inf_orig = Z_inf1_orig;
Malha(1).Z_med_orig = Z1_med_orig;

Malha(1).X = X1;
Malha(1).Y = Y1;
Malha(1).Z_sup = Z_sup1_orig; %Provavelmente não está correto, mas funciona para 
Malha(1).Z_inf = Z_inf1_orig; % a modelagem de algumas longarinas (H e I)
Malha(1).Z_med = Z1_med;

Malha(1).posMalha = posMalha;

%% INICIALIZANDO FIG %%
% figure()
% axis equal
% hold on
% axis tight
% grid on
% view(3)
% 
% % ================================ PLOT ================================= %
% surf(X1,Y1,Z_sup1,'FaceAlpha',0.3,'LineStyle',':')
% surf(X1,Y1,Z_inf1,'FaceAlpha',0.3)
% surf(X1,Y1,Z1_med)

end

%% Função aplicaRot() %%
% ============================== DESCRIÇÃO ============================== %
% Essa função recebe dois vetores referentes às coordenadas X e Z do perfil
% e aplica um twist.
function [ x_twist , z_sup_twist ] = aplicaRot( x , z_sup , twist , centro_twist )
x = x - centro_twist(1);
z_sup = z_sup - centro_twist(2);

R = [cosd(-twist) -sind(-twist); sind(-twist) cosd(-twist)];
v_sup = [x;z_sup];

v_sup_twist = R*v_sup;

x_twist = v_sup_twist(1,:);
z_sup_twist = v_sup_twist(2,:);

x_twist = x_twist + centro_twist(1);
z_sup_twist = z_sup_twist + centro_twist(2);

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

%% Função offsetMalha() %%
% ============================== DESCRIÇÃO ============================== %
% Essa função recebe uma malha de pontos e aplica um offset sobre ela.
function [ X_offset , Z_sup_offset , Z_inf_offset ] = offsetMalha( X , Z_sup , Z_inf , espessura )
[num_elementos_X,num_elementos_Y] = size(X);

X_offset = zeros(num_elementos_X,num_elementos_Y);
Z_sup_offset = X_offset;
Z_inf_offset = X_offset;

dx = diff(X);
dz_sup = diff(Z_sup);
dz_inf = diff(Z_inf);

vec_normal_sup = cat(3,dz_sup,-dx);
vec_normal_inf = cat(3,-dz_inf,dx);

versor_normal_sup = vec_normal_sup./sqrt(vec_normal_sup(:,:,1).^2+vec_normal_sup(:,:,2).^2);
versor_normal_inf = vec_normal_inf./sqrt(vec_normal_inf(:,:,1).^2+vec_normal_inf(:,:,2).^2);
versor_normal_sup = [versor_normal_sup;versor_normal_sup(end,:,:)];
versor_normal_inf = [versor_normal_inf;versor_normal_inf(end,:,:)];

for i = 1:num_elementos_Y    
    x_sup = X(:,i) + versor_normal_sup(:,i,1) * espessura;
    x_inf = X(:,i) + versor_normal_inf(:,i,1) * espessura;
    Z_sup_offset(:,i) = Z_sup(:,i) + versor_normal_sup(:,i,2) * espessura;
    Z_inf_offset(:,i) = Z_inf(:,i) + versor_normal_inf(:,i,2) * espessura;
    
    X_offset(:,i) = x_sup;
    Z_inf_offset(:,i) = interp1(x_inf,Z_inf_offset(:,i),x_sup,'linear','extrap');   
end
end

%% Função linspaceVet() %%
function y = linspaceVet(d1,d2,n)
n1 = n-1;
N = repmat(0:n1,length(d1),1);
y = d1 + N.*(d2-d1)./n1;
end

%% Função offsetBA_x() %%
function [ offset ] = offsetBA_x(Sb,l_BA)

offset = Sb .* tand(l_BA);
if length(offset) > 1
    for k = 2:1:length(offset)
        offset(k) = sum(offset(k-1:k));
    end
end
offset = [ 0 offset ];
end

%% Função lx() %%
function [ l_x ] = lx( Sb , c , l , n )
    l_x = zeros(1,length(Sb));
for i = 1:1:length(Sb)
    l_x(i) = (c(i)*(1/4-n(i)) + c(i+1)*(n(i+1)-1/4))/Sb(i) + tand(l(i));
    l_x(i) = atand(l_x(i));
end

end
%% Inverte Perfil (Para caso de caudas com perfis invertidos)
function[Z_sup1_orig,Z_inf1_orig] = ETT_InvertePerfil (Z_sup1_orig,Z_inf1_orig,X1_orig,Y1,x_BA,h_asa1_BA)
%Estabelecendo polo de rotacao
X_o = x_BA;
Y_o = 0;
Z_o = h_asa1_BA;
Origem_EH = [X_o;Y_o;Z_o]; %Origem do sistema de coordenadas da Empenagem e polo de rotacao

%Inializando variaveis auxiliares
Z_sup_linha = Z_sup1_orig;
X_sup_linha = X1_orig;
Y_sup_linha = Y1;

Z_inf_linha = Z_inf1_orig;
X_inf_linha = X1_orig;
Y_inf_linha = Y1;
%Passando as posicoes de uma base para outra
Z_sup = Z_sup1_orig - Z_o;
Z_inf = Z_inf1_orig - Z_o;
Y = Y1 - Y_o;
X = X1 - X_o;

%A inversão do perfil trata-se de uma rotação de 180º em relação ao eixo Y
%acompanhanda de uma rotação de 180º em relação ao eixo Z
ang_x = 0; %Angulo de rotação em X
ang_y = 180; %Angulo de rotação em Y
ang_z = 180; %Angulo de rotação em Z
rot_y = [cosd(ang_y) 0 -sind(ang_y); 0 1 0;sind(ang_y) 0 cosd(ang_y)];
rot_z = [cosd(ang_z) -sind(ang_z) 0; sind(ang_z) cosd(ang_z) 0; 0 0 1];

%Calculando
% [X_sup_linha;Y_sup_linha;Z_sup_linha] = rot_z.*rot_y.*[X;Y;Z_sup];
% 
% [X_inf_linha;Y_inf_linha;Z_inf_linha] = rot_z.*rot_y.*[X;Y;Z_inf];
    
%Retornando ao sistema de coordenadas global do avião



end