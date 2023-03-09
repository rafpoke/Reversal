%% Fun��o ETT_desRotElemento() %%
% ============================== Descri��o ============================== %
% Essa fun��o transforma os deslocamentos calculados em coordenadas globais
% para coordenadas locais.

function [ u_loc ] = ETT_desRotElemento( u , M_rot , elem )
%% Carregando Valores %%
conectividade = elem.conect;

num_elementos = size(conectividade,1);

%% Corpo da Fun��o %%
n_dim = 6;
n_nos = 2;
% ============================ Inicializa��o ============================ %
u_loc = zeros(n_nos*n_dim,1,num_elementos);
% ================== Transposi��o da matriz de rota��o ================== %
% Pra passar de coordenadas globais pra locais temos que utilizar a matriz
% de rota��o com �ngulos de Euler negativos. Isso equivale a transpor a
% matriz de rota��o.
% M_rot_transp = permute(M_rot,[2 1 3]);% t� errado mesmo isso aqui
% ================= Incrementa��o da matriz de rota��o ================== %
n_dim = 6;
n_nos = 2;
M_rot_inc = zeros(n_nos*n_dim,n_nos*n_dim,num_elementos);
M_rot_inc(1:3,1:3,:) = M_rot;
M_rot_inc(4:6,4:6,:) = M_rot;
M_rot_inc(7:9,7:9,:) = M_rot;
M_rot_inc(10:12,10:12,:) = M_rot;

% =============== Determina��o de u_loc pra cada elemento =============== %
for i = 1:num_elementos
    no1 = conectividade(i,1);
    no2 = conectividade(i,2);
    idx1 = [(no1-1)*n_dim+1:no1*n_dim];
    idx2 = [(no2-1)*n_dim+1:no2*n_dim];
    u_loc(:,1,i) = [ u(idx1) ; u(idx2) ];
    
    u_loc(:,1,i) = M_rot_inc(:,:,i) * u_loc(:,1,i);
end

end

