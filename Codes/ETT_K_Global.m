%% Função ETT_K_Global() %%
% ============================== Descrição ============================== %
% Essa função recebe as matrizes de rigidez locais (em coordenadas globais)
% e a matriz de conectividade e calcula a matriz de rigidez global do
% sistema.

function [K_glob] = ETT_K_Global( K_loc , nos , elem )
%% Carregando Valores %%
conectividade = elem.conect;
pontos = nos.pontos;

num_elem = size(conectividade,1);
num_nos = size(pontos,1);

%% Corpo da Função %%
n_dim = 6;
% ============================ Inicialização ============================ %
K_glob = zeros(num_nos*n_dim);
% K_glob = zeros(num_nos*n_dim,0) * zeros(0,num_nos*n_dim); que bizarro

% ============================= Scattering ============================== %
% Pra tentar entender isso aqui, dar uma olhada no Megson Cap. 6 (lá ele 
% faz na mão) e comparar. Buscar outras mídias também se possível, é chato 
% de entender. 

for i = 1:num_elem
    no1 = conectividade(i,1);
    no2 = conectividade(i,2);
    idx1 = (no1-1)*n_dim+1:no1*n_dim;
    idx2 = (no2-1)*n_dim+1:no2*n_dim;
    K_glob(idx1,idx1) = K_glob(idx1,idx1) + K_loc(1:n_dim,1:n_dim,i);
    K_glob(idx2,idx2) = K_glob(idx2,idx2) + K_loc(n_dim+1:2*n_dim,n_dim+1:2*n_dim,i);
    K_glob(idx1,idx2) = K_glob(idx1,idx2) + K_loc(1:n_dim,n_dim+1:2*n_dim,i);
    K_glob(idx2,idx1) = K_glob(idx1,idx2)';
end

% =================== Transformando em matriz esparsa =================== %
% Diminui drasticamente o tempo computacional
K_glob = sparse(K_glob);

end

