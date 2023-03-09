function [F_ext,struct_propriedades,struct_fem] = ETT_PreProcessingv7( aviao , P , F )
%% Carregando Valores %%
Longarina = aviao.Longarina;
Box = aviao.Box;
Tirante = aviao.Tirante;
% Cargas = aviao.Cargas;

D_ext_box = 4e-3;
D_int_box = 3e-3;
aviao.Iy_box = (pi/64) * (D_ext_box^4 - D_int_box^4);
aviao.Iz_box = (pi/64) * (D_ext_box^4 - D_int_box^4);
aviao.Ix_box = (pi/32) * (D_ext_box^4 - D_int_box^4);
aviao.A_box = (pi/4) * (D_ext_box^2 - D_int_box^2);

%% Definição dos nós %%
% Ordem dos nós na matriz (P são ops pontos de aplicação de força): 
% [ longarina inferior ; longarina superior ; longarina box ; P ];
% A ordenação dos nós do box é crescente no sentido vertical positivo

% % =========================== Forças aplicadas ========================== %
P = [ P{1} ; P{2} ; P{3} ];
F = [ F{1} ; F{2} ; F{3} ];
% ============================ Pontos do box ============================ %
pontos_box = Box.pontos;
% ========================= Pontos dos tirantes ========================= %
pontos_tirante1 = Tirante(1).pontos;
pontos_tirante2 = Tirante(2).pontos;
% ======================== Pontos das longarinas ======================== %
pontos_long = [ Longarina(1).pontos ; Longarina(2).pontos ];
idx_nos_long1 = 1:size(Longarina(1).pontos,1);
idx_nos_long2 = idx_nos_long1(end)+1:size(pontos_long,1);
idx_nos_long = 1:idx_nos_long2(end);
% ========================== Posição do engaste ========================= %
% engaste = [ 1 , 1 + size(Longarina(1).pontos,1) ];
engaste = 1;
% ===================== Pontos das longarinas + box ===================== %
% pontos_estrutura = unique([ pontos_long ; pontos_box ],'rows','stable');
% nos = [ pontos_long ; pontos_box(2:end-1,:) ];                 % faz a mesma coisa que a linha comentada só que mais rápido
nos = [ pontos_long ; pontos_box ];
idx_nos_box = idx_nos_long(end)+1:size(nos,1);
idx_nos_estrutura = 1:size(nos,1);
% =================== Pontos da estrutura + tirantes ==================== %
nos = [ nos ; pontos_tirante1 ; pontos_tirante2 ];
idx_nos_tirante1 = idx_nos_estrutura(end) + [ 1:size(pontos_tirante1,1) ];
idx_nos_tirante2 = idx_nos_tirante1(end) + [ 1:size(pontos_tirante2,1) ];
% ================= Nós da estrutura mais próximos de P ================= %
[X_E,X_P] = meshgrid(nos(idx_nos_long,1),P(:,1));
[Y_E,Y_P] = meshgrid(nos(idx_nos_long,2),P(:,2));
[Z_E,Z_P] = meshgrid(nos(idx_nos_long,3),P(:,3));

dist = sqrt((X_E-X_P).^2 + (Y_E-Y_P).^2 + (Z_E-Z_P).^2);
[~,idx_min_dist_P] = min(dist,[],2);
% ============= Nós da estrutura mais próximos dos tirantes ============= %
dist_inf = sqrt( (nos(idx_nos_long,1)-pontos_tirante1(1,1)).^2 + (nos(idx_nos_long,2)-pontos_tirante1(1,2)).^2 + (nos(idx_nos_long,3)-pontos_tirante1(1,3)).^2 );
dist_sup = sqrt( (nos(idx_nos_long,1)-pontos_tirante1(end,1)).^2 + (nos(idx_nos_long,2)-pontos_tirante1(end,2)).^2 + (nos(idx_nos_long,3)-pontos_tirante1(end,3)).^2 );
[~,idx_min_dist_inf] = min(dist_inf);
[~,idx_min_dist_sup] = min(dist_sup);
% ============================ Número de nós ============================ %
num_nos_box = length(idx_nos_box);
num_nos_tirante1 = length(idx_nos_tirante1);
num_nos_tirante2 = length(idx_nos_tirante2);
num_nos = size(nos,1);

%% Definição dos elementos %%
% ======================= Matriz de conectividade ======================= %
conectividade1 = [idx_nos_long1(1:end-1)',idx_nos_long1(2:end)'];
conectividade2 = [idx_nos_long2(1:end-1)',idx_nos_long2(2:end)'];
idx_elem_long1 = 1:size( conectividade1 , 1 );
idx_elem_long2 = size( conectividade1 , 1 )+1:size( conectividade1 , 1 )+size( conectividade2 , 1 );
idx_elem_long = 1:size( [ conectividade1 ; conectividade2 ] , 1 );

conectividade_box = flipud([ [idx_nos_long1(end),idx_nos_box(1)] ; [idx_nos_box(1:end-1)',idx_nos_box(2:end)'] ; [idx_nos_box(end),idx_nos_long2(end)]]);
idx_elem_box = idx_elem_long(end)+1:idx_elem_long(end)+size(conectividade_box,1);

conectividade_tirante1 = [ [idx_nos_estrutura(idx_min_dist_inf),idx_nos_tirante1(1)] ; [idx_nos_tirante1(1:end-1)',idx_nos_tirante1(2:end)'] ; [idx_nos_estrutura(idx_min_dist_sup),idx_nos_tirante1(end)] ];
conectividade_tirante2 = [ [idx_nos_estrutura(idx_min_dist_inf),idx_nos_tirante2(1)] ; [idx_nos_tirante2(1:end-1)',idx_nos_tirante2(2:end)'] ; [idx_nos_estrutura(idx_min_dist_sup),idx_nos_tirante2(end)] ];
idx_elem_tirante1 = size(  [ conectividade1 ; conectividade2 ; conectividade_box ] , 1 )+1:size(  [ conectividade1 ; conectividade2 ; conectividade_box ] , 1 )+size(conectividade_tirante1,1);
idx_elem_tirante2 = idx_elem_tirante1(end)+1:size([ conectividade1 ; conectividade2 ; conectividade_box ; conectividade_tirante1 ; conectividade_tirante2 ],1);
idx_elem_tirante = [ idx_elem_tirante1 , idx_elem_tirante2 ];

conectividade = [ conectividade1 ; conectividade2 ; conectividade_box ; conectividade_tirante1 ; conectividade_tirante2 ];
% ======================= Vetores de propriedades ======================= %
% Como as propriedades são dadas em cada nó no struct Longarina, mas
% queremos a propriedade em cada elemento da longarina (n nós --> n-1
% elementos) fazemos uma média entre os valores (resultando em um valor a
% menos que antes).
k = 10;
    % ================================ E ================================ %
    E1 = (Longarina(1).Ex(1:end-1)+Longarina(1).Ex(2:end))/2;
    E2 =(Longarina(2).Ex(1:end-1)+Longarina(2).Ex(2:end))/2;
    E_box = ones(1,num_nos_box+1) * Longarina(1).Ex(1);                     % CHUTE --> Box feito do mesmo material que a longarina
    E_t1 = [ k*E1(1) , (Tirante(1).Ex(1:end-1)+Tirante(1).Ex(2:end))/2 , k*E1(1) ];
    E_t2 = [ k*E1(1) , (Tirante(2).Ex(1:end-1)+Tirante(2).Ex(2:end))/2 , k*E1(1) ];
    E = [E1,E2,E_box,E_t1,E_t2];
    % ================================ G ================================ %
    G1 = (Longarina(1).Gx(1:end-1)+Longarina(1).Gx(2:end))/2;
    G2 = (Longarina(2).Gx(1:end-1)+Longarina(2).Gx(2:end))/2;
    G_box = ones(1,num_nos_box+1) * Longarina(1).Gx(1);                     % CHUTE --> Box feito do mesmo material que a longarina
    G_t1 = [ k*G1(1) , (Tirante(1).Gx(1:end-1)+Tirante(1).Gx(2:end))/2 , k*G1(1) ];
    G_t2 = [ k*G1(1) , (Tirante(2).Gx(1:end-1)+Tirante(2).Gx(2:end))/2 , k*G1(1) ];
    G = [G1,G2,G_box,G_t1,G_t2];
    % ================================ Ix =============================== %
    Ix1 = (Longarina(1).Iy(1:end-1)+Longarina(1).Iy(2:end))/2;
    Ix2 = (Longarina(2).Iy(1:end-1)+Longarina(2).Iy(2:end))/2;
    Ix_box = ones(1,num_nos_box+1) * aviao.Ix_box;
    Ix_t1 = [ k*Ix1(1) , (Tirante(1).Ix(1:end-1)+Tirante(1).Ix(2:end))/2 , k*Ix1(1) ];
    Ix_t2 = [ k*Ix1(1) , (Tirante(2).Ix(1:end-1)+Tirante(2).Ix(2:end))/2 , k*Ix1(1) ];
    Ix = [Ix1,Ix2,Ix_box,Ix_t1,Ix_t2];
    % ================================ Iy =============================== %
    Iy1 = (Longarina(1).Ix(1:end-1)+Longarina(1).Ix(2:end))/2;
    Iy2 = (Longarina(2).Ix(1:end-1)+Longarina(2).Ix(2:end))/2;
    Iy_box = ones(1,num_nos_box+1) * aviao.Iy_box;
    Iy_t1 = [ k*Iy1(1) , (Tirante(1).Iy(1:end-1)+Tirante(1).Iy(2:end))/2 , k*Iy1(1) ];
    Iy_t2 = [ k*Iy1(1) , (Tirante(2).Iy(1:end-1)+Tirante(2).Iy(2:end))/2 , k*Iy1(1) ];
    Iy = [Iy1,Iy2,Iy_box,Iy_t1,Iy_t2];
    % ================================ Iz =============================== %
    Iz1 = (Longarina(1).Iz(1:end-1)+Longarina(1).Iz(2:end))/2;
    Iz2 = (Longarina(2).Iz(1:end-1)+Longarina(2).Iz(2:end))/2;
    Iz_box = ones(1,num_nos_box+1) * aviao.Iz_box;
    Iz_t1 = [ k*Iz1(1) , (Tirante(1).Iz(1:end-1)+Tirante(1).Iz(2:end))/2 , k*Iz1(1) ];
    Iz_t2 = [ k*Iz1(1) , (Tirante(2).Iz(1:end-1)+Tirante(2).Iz(2:end))/2 , k*Iz1(1) ];
    Iz = [Iz1,Iz2,Iz_box,Iz_t1,Iz_t2];
    % ================================ A ================================ %
    A1 = (Longarina(1).A(1:end-1)+Longarina(1).A(2:end))/2;
    A2 = (Longarina(2).A(1:end-1)+Longarina(2).A(2:end))/2;
    A_box = ones(1,num_nos_box+1) * aviao.A_box;
    A_t1 = [ k*A1(1) , (Tirante(1).A(1:end-1)+Tirante(1).A(2:end))/2 , k*A1(1) ];
    A_t2 = [ k*A1(1) , (Tirante(2).A(1:end-1)+Tirante(2).A(2:end))/2 , k*A1(1) ];
    A = [A1,A2,A_box,A_t1,A_t2];
    % ================================ L ================================ %
    dist_nos = nos(conectividade(:,2),:) - nos(conectividade(:,1),:);
    L = sqrt(dist_nos(:,1).^2 + dist_nos(:,2).^2 + dist_nos(:,3).^2)';

%% Forças Externas %%
% ================================ F_ext ================================ %
% Ao invés de criar nós novos para os pontos de aplicação (altíssimo custo
% computacional), as forças aplicadas são transmitidas ao nó mais próximo
% da longarina (assim como um momento equivalente).
n_dim = 6;
F_ext = zeros(num_nos*n_dim,1);
Momento = cross( P - nos(idx_min_dist_P,:) , F );
for i = 1:size(P,1)
    F_ext( (idx_min_dist_P(i)-1)*n_dim + [1:n_dim] ) = F_ext( (idx_min_dist_P(i)-1)*n_dim + [1:n_dim] ) + [F(i,:)';Momento(i,:)'];
end

%% Simetria %%
% ================================= Nós ================================= %
nos_sim = nos;
nos_sim(:,2) = -nos_sim(:,2);
nos_sim(idx_nos_long2(1),:) = [];
% nos_sim([idx_nos_long1(1),idx_nos_long2(1)],:) = [];
nos = [ nos ; nos_sim ];
idx_nos_sim = num_nos+1:size(nos,1);
% ======================= Matriz de conectividade ======================= %
conectividade_sim = conectividade;
con_s = conectividade_sim;
con_s(con_s > idx_nos_long2(1)) = con_s(con_s > idx_nos_long2(1)) - 1;
con_s = con_s + num_nos;
con_s(idx_elem_long2(1),:) = [idx_nos_sim(idx_nos_long2(1)) , idx_nos_long2(1)];

conectividade = [ conectividade ; con_s ];
engaste = [ engaste , num_nos+1 ];
% num_nos_sim = size(nos_sim,1);
% conectividade_sim = [(conectividade1 + 1);(conectividade2 + 1);conectividade_box ; conectividade_tirante1 ; conectividade_tirante2 ] + num_nos_sim;
% conectividade_sim(idx_elem_long1(1),1) = idx_nos_long1(1);
% conectividade_sim(idx_elem_long2(1),1) = idx_nos_long2(1);
% conectividade = [ conectividade ; conectividade_sim ];
% ======================= Vetores de propriedades ======================= %
E = [E,E];
G = [G,G];
Ix = [Ix,Ix];
Iy = [Iy,Iy];
Iz = [Iz,Iz];
A = [A,A];
L = [L,L];
% ================================ F_ext ================================ %
F_ext_sim = F_ext;
vet_nos = 1:num_nos;
idx_F = [(vet_nos-1)*n_dim+2,(vet_nos-1)*n_dim+4,(vet_nos-1)*n_dim+6];
F_ext_sim(idx_F) = - F_ext_sim(idx_F);
F_ext_sim((idx_nos_long2(1)-1)*n_dim+(1:6),:) = [];
F_ext = [ F_ext ; F_ext_sim ];

% elem = num2cell((1:size(conectividade,1))');
% table(elem,conectividade)

%     elem = num2cell([1:size(conectividade,1)]');
%     table(elem,E',G',Ix',Iy',Iz',A',L')

%% Outputs %%
% ========================= struct_propriedades ========================= %
struct_propriedades = struct;
struct_propriedades.E = E;
struct_propriedades.G = G;
struct_propriedades.Ix = Ix;
struct_propriedades.Iy = Iy;
struct_propriedades.Iz = Iz;
struct_propriedades.A = A;
struct_propriedades.L = L;
% ============================= struct_fem ============================== %
struct_fem = struct;
struct_fem.nos = nos;
struct_fem.conectividade = conectividade;
struct_fem.engaste = engaste;
struct_fem.idx_nos_long1 = idx_nos_long1;
struct_fem.idx_nos_long2 = idx_nos_long2;
struct_fem.idx_nos_box = idx_nos_box;
struct_fem.idx_nos_estrutura = idx_nos_estrutura;
struct_fem.idx_nos_tirante1 = idx_nos_tirante1;
struct_fem.idx_nos_tirante2 = idx_nos_tirante2;
struct_fem.idx_elem_long = idx_elem_long;
struct_fem.idx_elem_long1 = idx_elem_long1;
struct_fem.idx_elem_long2 = idx_elem_long2;
struct_fem.idx_elem_box = idx_elem_box;
struct_fem.idx_elem_tirante1 = idx_elem_tirante1;
struct_fem.idx_elem_tirante2 = idx_elem_tirante2;
struct_fem.idx_elem_tirante = idx_elem_tirante;

end

%% Função linspaceVet() %%
function y = linspaceVet(d1,d2,n)
n1 = n-1;
N = repmat(0:n1,length(d1),1);
y = d1 + N.*(d2-d1)./n1;
end
