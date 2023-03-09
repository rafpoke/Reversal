%% Função ETT_Solve_KuF() %%
% ============================== Descrição ============================== %
% Essa função resolve o sistema F = K*u. São aplicadas também as condições
% de contorno u(cc) = 0. Geralmente em textos de MEF as condições de
% contorno são aplicadas através de uma matriz que você tem que subtrair.
% Por algum motivo o jeito que eu faço de resolver o sistema sem as linhas
% e colunas que estão sendo restringidas funciona. Eu não sei o porquê,
% descobri isso meio na sorte e funcionou.

function [u] = ETT_Solve_KuF( K_glob , F_ext , nos )
%% Carregando Valores %%
cc = nos.cc;
pontos = nos.pontos;

num_nos = size(pontos,1);

%% Corpo da Função %%
n_dim = 6;
% =============================== F = Ku ================================ %
u = zeros(num_nos*n_dim,1);
u(~cc) = K_glob(~cc,~cc)\F_ext(~cc);

end

