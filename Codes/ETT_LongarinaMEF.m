%% Fun��o ETT_LongarinaMEF() %%
% ============================== Descri��o ============================== %
% Essa fun��o implementa uma an�lise de em MEF da longarina. Todas as
% estruturas s�o aproximados por elementos de barra discretos. 

% Inicialmente � realizado um pr�-processamento, que consiste em descrever 
% a estrutura a partir de n�s e elementos de barra e calcular as for�as e 
% momentos agindo sobre cada n�. Em seguida � calculada as matrizes de 
% cossenos diretores de cada elemento da malha que permite transformar as 
% matrizes de rigidez de um referencial local para o referencial global.

% O processamento consiste em calcular as matrizes de rigidez locais de
% cada elemento e realizar a transforma��o para uma matriz de rigidez
% global do sistema. Em seguida resolve-se o sistema linear F = K*u,
% descobrindo os deslocamentos do n�s (em coordenadas globais).

% O p�s-processamento transforma esses deslocamentos globais em
% deslocamentos locais (no sistema de coordenadas de cada elemento). Dessa
% forma pode-se aplicar as equa��es (6.44) do Megson para determinar as
% for�as cortantes e de tra��o assim como momentos fletores e torsores
% atuantes em cada elemento. 

% =============================== Outputs =============================== %
% - Longarina - struct [kx1] (onde k � o n�mero de "longarinas" da
% aeronave. Incluem-se nesse lista quaisquer outras estruturas adicionais
% da aeronave relevantes para a an�lise, como tirantes, stringers, nervuras
% etc). Essa struct cont�m todas as informa��es relevantes de cada
% longarina, como: matriz de pontos, vetor de in�rcias, informa��es
% geom�tricas, massa etc.

% - F_int - matriz [elx3] (onde el � o n�mero de elementos da malha). A
% i-�sima linha da matriz descreve as for�as internas [Fx,Fy,Fz] 
% (em coordenadas locais) no i-�simo elemento (o i-�simo elementos � aquele
% descrito pela i-�sima linha da matriz de conectividade elem.conect).

% - M_int - matriz [elx3] (onde el � o n�mero de elementos da malha). A
% i-�sima linha da matriz descreve os momentos internos [Mx,My,Mz] 
% (em coordenadas locais) no i-�simo elemento (o i-�simo elementos � aquele
% descrito pela i-�sima linha da matriz de conectividade elem.conect).

% - u - vetor [n_dim*nx1] (onde n � o n�mero de n�s da malha e n_dim � o
% n_dim � o n�mero de graus de liberdade de cada n�, nesse caso n_dim = 6).
% Esse vetor indica os deslocamentos (em coordenadas globais) dos n�s da
% malha. Os deslocamentos s�o ordenados da seguinte forma: 
% [dx;dy;dz;tx;ty;tz] onde d indica deslocamento linear e t deslocamento
% angular.

% - nos - struct que descreve os n�s da malha. Cont�m os seguintes campos:
%   - pontos - matriz [nx3] (onde n � o n�mero de n�s). Cont�m as
%   coordenadas [x,y,z] de cada n�.
%   - idx - vetor [1xn] (onde n � o n�mero de n�s). Relaciona cada n� da
%   estrutura com uma longarina. 
%   Ex: idx = [1 1 1 2 2] significa que os n�s [1,2,3] pertencem �
%   longarina 1, enquanto que os n�s [4,5] pertencem � longarina 2.
%   OBS IMPORTANTE: assume-se no c�digo que a sequ�ncia de n�meros ser�
%   sempre crescente e sem repeti��es (os n�s de cada longarina est�o
%   agrupados juntos). Ou seja, algo como idx = [1 1 1 2 2 1] nunca pode
%   acontecer.
%   - cc - vetor l�gico [n_dim*nx1] (onde n � o n�mero de n�s da malha e 
%   n_dim � o n�mero de graus de liberdade de cada n�, nesse caso n_dim =
%   6). Esse � um vetor de condi��es de contorno que indica quais
%   deslocamentos do vetor u s�o restringidos.
% Nota-se que a numera��o dos n�s � dada pela matriz de pontos e �
% respeitada nos outros campos. Ou seja, a i-�sima linha da matriz de
% pontos se refere aos pontos do i-�simo n�, o i-�simo termo de idx
% refere-se � longarina da qual o i-�simo n� faz parte e as linhas
% ((i-1)*n_dim + [1:6]) de cc referem-se �s condi��es de contorno
% [dx,dy,dz,tx,ty,tz] do i-�simo n�.

% - elem - struct que descreve os elementos da malha. Cont�m os seguintes
% campos:
%   - conect - matriz [elx2] (onde el � o n�mero de elementos da malha).
%   � a chamada matriz de conectividade que define como que os n�s se
%   conectam para formar os elementos. Cada linha indica os dois n�s que
%   formam um elemento. Ex: [1 2;2 3;3 1] descreve uma malha com 3
%   elementos que formam um tri�ngulo (n� 1 conecta em 2, que conecta em 3,
%   que conecta em 1 de novo).
%   - idx - vetor [elx1] (onde el � o n�mero de elementos). Relaciona cada
%   elemento com uma longarina.
%   Ex: idx = [1 1 1 2 2] significa que os elementos [1,2,3] pertencem �
%   longarina 1, enquanto que os elementos [4,5] pertencem � longarina 2.
%   - A,Ix,Iy,Iz,G,E,L vetores [1xel] (onde el � o n�mero de elementos da
%   malha). Esses vetores descrevem, respectivamente, �rea, momento de
%   in�rcia em x (coordenadas locais), momento de in�rcia em y (coordenadas
%   locais), momentos de in�rcia em z (coordenadas locais), m�dulo de
%   cisalhamento, m�dulo de elasticidade e comprimento de cada elemento.
% Nota-se que a numera��o dos n�s � dada pela matriz de conectividade e �
% respeitada pelos outros campos. Ou seja, a i-�sima linha da matriz de
% conectividade indica os n�s que comp�em o i-�simo elemento, o i-�simo
% termo de ids refere-se � longarina da qual o i-�simo elemento faz parte e
% o i-�simo termo de (A,Ix,Iy,Iz,G,E,L) refere-se � respectiva propriedade
% do i-�simo elemento.


% A fun��o d� erro se a barra estiver alinhada com o eixo Z (vetor de
% deslocamento NaN).

function [F_int,M_int,u,nos,elem,Longarina] = ETT_LongarinaMEF( Longarina , F , P , varargin )
%% Pre-Processing %%
% ================== Vetor de for�as externas e structs ================= %
[F_ext,nos,elem] = ETT_PreProcessing( Longarina , F , P );
% ==================== Matrizes de cossenos diretores =================== %
[M_rot,elem] = ETT_CossenosDiretores( nos , elem );

%% Processing %%
% =========================== Matriz K local ============================ %
K_loc = ETT_K_Local( elem );
% ========================== Rota��o de K_loc =========================== %
K_loc_rot = ETT_Rot_K_Local( K_loc , M_rot );
% =========================== Matriz K Global =========================== %
K_glob = ETT_K_Global( K_loc_rot , nos , elem );
% ================== Resolu��o do sistema linear Ku=F =================== %
u = ETT_Solve_KuF( K_glob , F_ext , nos );

%% Post-Processing %%
% ========================== Rea��es Internas =========================== %
u_loc = ETT_desRotElemento( u , M_rot , elem );
[F_int,M_int] = ETT_ForcasIntElemento( u_loc , elem );
% ============================= Frescurinhas ============================ %
% tableU(u,F_int,M_int,F_ext,elem);
% tableElem( elem );
% % mostarK(K_loc_rot,K_glob,struct_fem)
% plotDesloc(elem,nos,u,F_ext)
Longarina = ETT_F_FEM( Longarina , nos , F_ext );

end

