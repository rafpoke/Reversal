%% Fun��o ETT_Solve_KuF() %%
% ============================== Descri��o ============================== %
% Essa fun��o resolve o sistema F = K*u. S�o aplicadas tamb�m as condi��es
% de contorno u(cc) = 0. Geralmente em textos de MEF as condi��es de
% contorno s�o aplicadas atrav�s de uma matriz que voc� tem que subtrair.
% Por algum motivo o jeito que eu fa�o de resolver o sistema sem as linhas
% e colunas que est�o sendo restringidas funciona. Eu n�o sei o porqu�,
% descobri isso meio na sorte e funcionou.

function [u] = ETT_Solve_KuF( K_glob , F_ext , nos )
%% Carregando Valores %%
cc = nos.cc;
pontos = nos.pontos;

num_nos = size(pontos,1);

%% Corpo da Fun��o %%
n_dim = 6;
% =============================== F = Ku ================================ %
u = zeros(num_nos*n_dim,1);
u(~cc) = K_glob(~cc,~cc)\F_ext(~cc);

end

