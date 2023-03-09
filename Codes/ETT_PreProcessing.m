%% Função ETT_PreProcessing() %%
% ============================== Descrição ============================== %
% Essa função é responsável por receber a struct de longarina, as forças
% aplicadas e seus pontos de aplicação e gerar as estruturas de dados
% necessárias para realizar as computações de elementos finitos. Isso
% inclui:
% - Criação de nós - descrever o pontos de cada nó e a qual longarina eles
% pertencem.
% - Criação de elementos - descrever as longarinas a partir de elementos
% (união de dois nós). Descrever os nós que compõem cada elemento, assim
% como suas propriedades e a qual longarina esses elementos pertencem.
% - Condições de contorno - definir quais graus de liberdade de quais nós
% são restringidos a partir de um vetor.
% - Forças externas - transladar as forças aplicadas nos pontos P para os
% nós e aplicar momentos equivalentes que compensem essa translação.

function [F_ext,nos,elem] = ETT_PreProcessing( Longarina , F , P )
%% Carregando Valores %%
L = Longarina;
num_pontos = [L.num_pontos];
cc = [L.cc];
connect = {L.connect};

%% Definição do Nós %%
nos = struct;
% ========================== Pontos Longarina =========================== %
Pt = arrayfun(@(x) x.pontos',L,'UniformOutput',false);
nos.pontos = cell2mat(Pt)';

% ============================== idx_long =============================== %
nos.idx = spread(cumsum(num_pontos));                                       % a qual longarina pertence cada nó

% % ========================== Forças Aplicadas =========================== %
% F_ext = fext( nos , L , P , F );                                            % forças e momentos aplicados a cada nó

%% Condições de Contorno %%
% =========================== Separando Cell ============================ %
idx_long_cc = [cc{1,:}];                                                    % indexamento de quais longarina que possuem condições de contorno
idx_nos_loc = [cc{2,:}];                                                    % indexamento local (internamente em cada longarina) do nó que está sendo restringido
idx_cc = {cc{3,:}};                                                         % indexamento dos graus de liberdade restringidos (x,y,z,tx,ty,tz)

% ==================== Transformação pra Nós Globais ==================== %
aux = [0 cumsum([L.num_pontos])];
idx_nos_glob = idx_nos_loc + aux( idx_long_cc );

% ============ Indexamento dos Graus de Liberdade Restritos ============= %
n_dim = 6;
idx_gl = zeros( n_dim * size(nos.pontos,1) , 1 );

for i = 1:length(idx_nos_glob)
    idx_gl( (idx_nos_glob(i)-1)*n_dim + idx_cc{i} ) = 1;
end

nos.cc = logical( idx_gl );                                                 % indexamento global dos graus de liberdade restritos

%% Definição dos Elementos %%
elem = struct;
% ==================== Conectividade Intra Estrutura ==================== %
% Aqui são definidos os elementos de cada longarina separadamente.
[elem.conect , elem.idx] = conect_intra( nos );

% ==================== Conectividade Inter Estrutura ==================== %
% Aqui são definidos os elementos que conectam as longarinas entre si.
nos.conect = connect;
[elem.conect , elem.idx , L ] = conect_inter( nos , elem , L );

% ============================= Propriedades ============================ %
elem = elem_prop( elem , nos , L );

%% Tira Nós Iguais e Elementos Iguais %%

pontos_round = round( nos.pontos , 8 );
% [nos.pontos,IA,IC] = uniquetol(nos.pontos,'ByRows',true);
[~,IA,IC] = unique(pontos_round,'rows','stable');
nos.pontos = nos.pontos(IA,:);
[nos,elem] = atualizaStructs(nos,elem,IA,IC);

%% Forças %%
% ========================== Forças Aplicadas =========================== %
F_ext = fext( nos , L , P , F );                                            % forças e momentos aplicados a cada nó

%% PLOT %%

% figure(1)
% hold on
% axis equal
% 
% plotF( nos , F_ext );
% plotElem( nos , elem );

end

%% Função atualizaStructs() %%
function [ nos , elem ] = atualizaStructs(nos,elem,idx_new,idx_old)

% ============================== Nos novo =============================== %
n_dim = 6;
% idx_gl = 1:(length(idx_new))*n_dim;
idx_gl = ((idx_new-1)*n_dim + (1:6))';
idx_gl = idx_gl(:);

nos.idx = nos.idx( idx_new );
nos.cc = nos.cc( idx_gl );

% ========================== Conectividade nova ========================= %
conect_new = zeros(size(elem.conect));

for i = 1:length(idx_old)
    conect_new( elem.conect(:) == i ) = idx_old(i); 
end

elem.conect = conect_new;

% ================ Retirada de elementos de tamanhho nulo =============== %
err = 1e-10;

elem_null = elem.L < err;

elem.conect(elem_null,:) = [];
elem.idx(elem_null) = [];
elem.A(elem_null) = [];
elem.Ix(elem_null) = [];
elem.Iy(elem_null) = [];
elem.Iz(elem_null) = [];
elem.E(elem_null) = [];
elem.G(elem_null) = [];
elem.L(elem_null) = [];

end

%% Função fext() %%
function [ F_ext ] = fext( nos , L , P , F )
% ========================= Longarinas Elegíveis ======================== %
idx_F = find( [L.idx_F] == true );

% ============================ Nós Elegíveis ============================ %
idx_nos_ele = logical( sum( nos.idx == idx_F' , 1 ) );
nos.pontos(~idx_nos_ele,:) = nan;

% ======================= Nós Mais Próximos de P ======================== %
% [X_E,X_P] = meshgrid(nos.pontos(idx_nos_ele,1),P(:,1));
% [Y_E,Y_P] = meshgrid(nos.pontos(idx_nos_ele,2),P(:,2));
% [Z_E,Z_P] = meshgrid(nos.pontos(idx_nos_ele,3),P(:,3));

[X_E,X_P] = meshgrid(nos.pontos(:,1),P(:,1));
[Y_E,Y_P] = meshgrid(nos.pontos(:,2),P(:,2));
[Z_E,Z_P] = meshgrid(nos.pontos(:,3),P(:,3));

dist = sqrt((X_E-X_P).^2 + (Y_E-Y_P).^2 + (Z_E-Z_P).^2);
% dist = (Y_E-Y_P).^2 + (Z_E-Z_P).^2;
[~,idx_min_dist_P] = min(dist,[],2);

% k = 10;
% [~,idx_min_dist_P] = sort(dist,2);
% idx_min_dist_P = idx_min_dist_P(:,1:k);

% % ================================ F_ext ================================ %
% num_nos = size(nos.pontos,1);
% n_dim = 6;
% F_ext = zeros(num_nos*n_dim,1);
% for j = 1:k
%     Momento = cross( P - nos.pontos(idx_min_dist_P(:,j),:) , F );
%     for i = 1:size(P,1)
%         F_ext( (idx_min_dist_P(i,j)-1)*n_dim + [1:n_dim] ) = F_ext( (idx_min_dist_P(i,j)-1)*n_dim + [1:n_dim] ) + [F(i,:)';Momento(i,:)']/k;
%     end
% end
% 
% sum(F_ext)

num_nos = size(nos.pontos,1);
n_dim = 6;
F_ext = zeros(num_nos*n_dim,1);
Momento = cross( P - nos.pontos(idx_min_dist_P,:) , F );
for i = 1:size(P,1)
    F_ext( (idx_min_dist_P(i)-1)*n_dim + [1:n_dim] ) = F_ext( (idx_min_dist_P(i)-1)*n_dim + [1:n_dim] ) + [F(i,:)';Momento(i,:)'];
end

end

%% Função elem_prop() %%
function [ elem ] = elem_prop( elem , nos , L )

dist = nos.pontos(elem.conect(:,2),:) - nos.pontos(elem.conect(:,1),:);

% Isso aqui assume que as longarinas vão estar agrupadas, ou seja, não vai
% ter metade de uma longarina em um lugar da matriz de conectividade e
% outra metade em outro lugar. Pode ser interessante criar um erro pra caso
% essa condição seja quebrada.
aux = [elem.idx(diff(elem.idx) ~= 0)' elem.idx(end)];

elem.A = [L(aux).A];

elem.Ix = [L(aux).Ix];
elem.Iy = [L(aux).Iy];
elem.Iz = [L(aux).Iz];

% elem.Ex = [L(aux).Ex];
% elem.Ey = [L(aux).Ey];
% elem.Ez = [L(aux).Ez];
elem.E = [L(aux).Ex];

% elem.Gx = [L(aux).Gx];
% elem.Gy = [L(aux).Gy];
% elem.Gz = [L(aux).Gz];
elem.G = [L(aux).Gx];

elem.L = sqrt( sum( dist.^2 , 2 ) )';

end

%% Função conect_intra() %%
function [ conect , conect_map ] = conect_intra( nos )

idx_nos = nos.idx;

t = (1:length(idx_nos))';
conect = [t(1:end-1),t(2:end)];

conect_map = idx_nos(conect);
conect( conect_map(:,1) ~= conect_map(:,2) , : ) = [];
conect_map( conect_map(:,1) ~= conect_map(:,2) , : ) = [];
conect_map(:,2) = [];

end

%% Função conect_inter() %%
function [ conect , conect_map , L ] = conect_inter( nos , elem , L )

conect = elem.conect;
conect_map = elem.idx; 

cin = nos.conect;
idx_nos = nos.idx;

for i = 1:length( cin )
    if ~ isempty( cin{i} )
        for j = 1:size( cin{i} , 1 )
        
            long1 = cin{i}{j,1};
            long2 = cin{i}{j,2};    
            method_no1 = cin{i}{j,3};
            method_no2 = cin{i}{j,4};

            is_long1 = idx_nos == long1;
            is_long2 = idx_nos == long2;

            no1_loc = method_no1( L(long1) );
            no1_glob = loc2glob( no1_loc , idx_nos , long1 );

            if method_no2 == 'nearest'
                nos_long2 = nos.pontos( is_long2 , : );
                ponto_no1 = nos.pontos( no1_glob , : );

                dist = ponto_no1 - nos.pontos( is_long2 , : );
                [~,no2_loc] = min( sum( dist.^2 , 2 ) );

                no2_glob = loc2glob( no2_loc , idx_nos , long2 );
            end 

            conect = [conect; no1_glob no2_glob];
            conect_map(end+1) = max(conect_map) + 1;

            L = defaultElementosInter( L );
        
        end
    end   
end

end

%% Função defaultElementosInter() %%
function [ L ] = defaultElementosInter( L )

idx = length(L) + 1;
% L(idx).A = 25e-6;
% L(idx).Ix = 2e-8;
% L(idx).Iy = 2e-9;
% L(idx).Iz = 2e-9;
% L(idx).Ex = 7e10;
% L(idx).Gx = 5e9;

L(idx).A = 25e-6;
L(idx).Ix = 20e-8;
L(idx).Iy = 20e-9;
L(idx).Iz = 20e-9;
L(idx).Ex = 7e10;
L(idx).Gx = 5e9;

end

%% Função loc2glob() %%
function [idx_glob] = loc2glob( idx_loc , idx_long , num_long )

first_no_long = find( idx_long == num_long , 1 , 'first' );
idx_glob = first_no_long + idx_loc - 1;

end

%% Função spread() %%
function [ new_array ] = spread( points )

points = [0 points];

new_array = zeros(1,points(end));

for i = 1:length(points)-1
    new_array(points(i)+1:points(i+1)) = i; 
end

end

%% Função plotF %%
function [ ] = plotF( nos , F_ext )

P = nos.pontos;

n_dim = 6;
num_nos = length(F_ext)/n_dim;

aux = 0:num_nos-1;

F = F_ext( aux * n_dim + [1;2;3] )'; F = F / max(abs(F(:)));
M = F_ext( aux * n_dim + [4;5;6] )'; M = M / max(abs(M(:)));

for i = 1:size(P,1)
    qF = quiver3(P(i,1),P(i,2),P(i,3),F(i,1),F(i,2),F(i,3)); qF.Color = [235, 59, 90]/256;
%     qM = quiver3(P(i,1),P(i,2),P(i,3),M(i,1),M(i,2),M(i,3)); qM.Color = [165, 94, 234]/256;
end

end

%% Função plotElem() %%
function [] = plotElem( nos , elem )

conectividade = elem.conect;

for k = 1:size(conectividade,1)
%     i = idx(k);
    i = k;
    inicio = nos.pontos(conectividade(i,1),:);
    fim = nos.pontos(conectividade(i,2),:);
    linha = [inicio;fim];
    plot3(linha(:,1),linha(:,2),linha(:,3),'LineWidth',2,'Color',[0.4 0.4 0.4])
    hold on
    axis equal
    grid on
end

end