%% Função ETT_LongarinaMEF() %%
% ============================== Descrição ============================== %
% Essa função implementa uma análise de em MEF da longarina. Todas as
% estruturas são aproximados por elementos de barra discretos. 

% Inicialmente é realizado um pré-processamento, que consiste em descrever 
% a estrutura a partir de nós e elementos de barra e calcular as forças e 
% momentos agindo sobre cada nó. Em seguida é calculada as matrizes de 
% cossenos diretores de cada elemento da malha que permite transformar as 
% matrizes de rigidez de um referencial local para o referencial global.

% O processamento consiste em calcular as matrizes de rigidez locais de
% cada elemento e realizar a transformação para uma matriz de rigidez
% global do sistema. Em seguida resolve-se o sistema linear F = K*u,
% descobrindo os deslocamentos do nós (em coordenadas globais).

% O pós-processamento transforma esses deslocamentos globais em
% deslocamentos locais (no sistema de coordenadas de cada elemento). Dessa
% forma pode-se aplicar as equações (6.44) do Megson para determinar as
% forças cortantes e de tração assim como momentos fletores e torsores
% atuantes em cada elemento. 

% =============================== Outputs =============================== %
% - Longarina - struct [kx1] (onde k é o número de "longarinas" da
% aeronave. Incluem-se nesse lista quaisquer outras estruturas adicionais
% da aeronave relevantes para a análise, como tirantes, stringers, nervuras
% etc). Essa struct contém todas as informações relevantes de cada
% longarina, como: matriz de pontos, vetor de inércias, informações
% geométricas, massa etc.

% - F_int - matriz [elx3] (onde el é o número de elementos da malha). A
% i-ésima linha da matriz descreve as forças internas [Fx,Fy,Fz] 
% (em coordenadas locais) no i-ésimo elemento (o i-ésimo elementos é aquele
% descrito pela i-ésima linha da matriz de conectividade elem.conect).

% - M_int - matriz [elx3] (onde el é o número de elementos da malha). A
% i-ésima linha da matriz descreve os momentos internos [Mx,My,Mz] 
% (em coordenadas locais) no i-ésimo elemento (o i-ésimo elementos é aquele
% descrito pela i-ésima linha da matriz de conectividade elem.conect).

% - u - vetor [n_dim*nx1] (onde n é o número de nós da malha e n_dim é o
% n_dim é o número de graus de liberdade de cada nó, nesse caso n_dim = 6).
% Esse vetor indica os deslocamentos (em coordenadas globais) dos nós da
% malha. Os deslocamentos são ordenados da seguinte forma: 
% [dx;dy;dz;tx;ty;tz] onde d indica deslocamento linear e t deslocamento
% angular.

% - nos - struct que descreve os nós da malha. Contém os seguintes campos:
%   - pontos - matriz [nx3] (onde n é o número de nós). Contém as
%   coordenadas [x,y,z] de cada nó.
%   - idx - vetor [1xn] (onde n é o número de nós). Relaciona cada nó da
%   estrutura com uma longarina. 
%   Ex: idx = [1 1 1 2 2] significa que os nós [1,2,3] pertencem à
%   longarina 1, enquanto que os nós [4,5] pertencem à longarina 2.
%   OBS IMPORTANTE: assume-se no código que a sequência de números será
%   sempre crescente e sem repetições (os nós de cada longarina estão
%   agrupados juntos). Ou seja, algo como idx = [1 1 1 2 2 1] nunca pode
%   acontecer.
%   - cc - vetor lógico [n_dim*nx1] (onde n é o número de nós da malha e 
%   n_dim é o número de graus de liberdade de cada nó, nesse caso n_dim =
%   6). Esse é um vetor de condições de contorno que indica quais
%   deslocamentos do vetor u são restringidos.
% Nota-se que a numeração dos nós é dada pela matriz de pontos e é
% respeitada nos outros campos. Ou seja, a i-ésima linha da matriz de
% pontos se refere aos pontos do i-ésimo nó, o i-ésimo termo de idx
% refere-se à longarina da qual o i-ésimo nó faz parte e as linhas
% ((i-1)*n_dim + [1:6]) de cc referem-se às condições de contorno
% [dx,dy,dz,tx,ty,tz] do i-ésimo nó.

% - elem - struct que descreve os elementos da malha. Contém os seguintes
% campos:
%   - conect - matriz [elx2] (onde el é o número de elementos da malha).
%   É a chamada matriz de conectividade que define como que os nós se
%   conectam para formar os elementos. Cada linha indica os dois nós que
%   formam um elemento. Ex: [1 2;2 3;3 1] descreve uma malha com 3
%   elementos que formam um triângulo (nó 1 conecta em 2, que conecta em 3,
%   que conecta em 1 de novo).
%   - idx - vetor [elx1] (onde el é o número de elementos). Relaciona cada
%   elemento com uma longarina.
%   Ex: idx = [1 1 1 2 2] significa que os elementos [1,2,3] pertencem à
%   longarina 1, enquanto que os elementos [4,5] pertencem à longarina 2.
%   - A,Ix,Iy,Iz,G,E,L vetores [1xel] (onde el é o número de elementos da
%   malha). Esses vetores descrevem, respectivamente, área, momento de
%   inércia em x (coordenadas locais), momento de inércia em y (coordenadas
%   locais), momentos de inércia em z (coordenadas locais), módulo de
%   cisalhamento, módulo de elasticidade e comprimento de cada elemento.
% Nota-se que a numeração dos nós é dada pela matriz de conectividade e é
% respeitada pelos outros campos. Ou seja, a i-ésima linha da matriz de
% conectividade indica os nós que compõem o i-ésimo elemento, o i-ésimo
% termo de ids refere-se à longarina da qual o i-ésimo elemento faz parte e
% o i-ésimo termo de (A,Ix,Iy,Iz,G,E,L) refere-se à respectiva propriedade
% do i-ésimo elemento.


% A função dá erro se a barra estiver alinhada com o eixo Z (vetor de
% deslocamento NaN).

function [F_int,M_int,u,nos,elem,Longarina] = ETT_LongarinaMEF( Longarina , F , P , varargin )
%% Pre-Processing %%
% ================== Vetor de forças externas e structs ================= %
[F_ext,nos,elem] = ETT_PreProcessing( Longarina , F , P );
% ==================== Matrizes de cossenos diretores =================== %
[M_rot,elem] = ETT_CossenosDiretores( nos , elem );

%% Processing %%
% =========================== Matriz K local ============================ %
K_loc = ETT_K_Local( elem );
% ========================== Rotação de K_loc =========================== %
K_loc_rot = ETT_Rot_K_Local( K_loc , M_rot );
% =========================== Matriz K Global =========================== %
K_glob = ETT_K_Global( K_loc_rot , nos , elem );
% ================== Resolução do sistema linear Ku=F =================== %
u = ETT_Solve_KuF( K_glob , F_ext , nos );

%% Post-Processing %%
% ========================== Reações Internas =========================== %
u_loc = ETT_desRotElemento( u , M_rot , elem );
[F_int,M_int] = ETT_ForcasIntElemento( u_loc , elem );
% ============================= Frescurinhas ============================ %
% tableU(u,F_int,M_int,F_ext,elem);
% tableElem( elem );
% % mostarK(K_loc_rot,K_glob,struct_fem)
% plotDesloc(elem,nos,u,F_ext)
Longarina = ETT_F_FEM( Longarina , nos , F_ext );

end

