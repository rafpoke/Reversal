%% Função ETT_Rot_K_Local() %%
% ============================== Descrição ============================== %
% Essa função recebe as matrizes de rigidez locais (em coordenadas locais)
% e as matrizes de rotação de cada elemento e calcula as matrizes de
% rigidez locais em coordenadas globais.

function [K_loc] = ETT_Rot_K_Local( K_loc , M_rot )
%% Carregando Valores %%
num_elementos = size(M_rot,3);

%% Corpo da Função %%
% ================= Incrementação da matriz de rotação ================== %
n_dim = 6;
n_nos = 2; %número de nós em cada elemento 
M_rot_inc = zeros(n_nos*n_dim,n_nos*n_dim,num_elementos);
M_rot_inc(1:3,1:3,:) = M_rot;
M_rot_inc(4:6,4:6,:) = M_rot;
M_rot_inc(7:9,7:9,:) = M_rot;
M_rot_inc(10:12,10:12,:) = M_rot;
% ================== Transposição da matriz de rotação ================== %
M_rot_inc_T = permute(M_rot_inc,[2 1 3]);
% =================== Cálculo da matriz K rotacionada =================== %
for i = 1:num_elementos
    K_loc(:,:,i) = M_rot_inc_T(:,:,i) * K_loc(:,:,i) * M_rot_inc(:,:,i);
end

end

