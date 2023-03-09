%% Função ETT_K_Local() %%
% ============================== Descrição ============================== %
% Essa função calcula a matriz de rigidez local de cada elemento (em
% coordenada locais). O formato de matriz utilizado é o de um elemento de
% barra com 12 graus de liberdade.

function [K] = ETT_K_Local( elem )
% ========================== Carregando struct ========================== %
E = elem.E;
G = elem.G;
Ix = elem.Ix;
Iy = elem.Iy;
Iz = elem.Iz;
A = elem.A;
L = elem.L;
% ======================= Transposição dos vetores ====================== %
E = reshape(E,1,1,[]);
G = reshape(G,1,1,[]);
Ix = reshape(Ix,1,1,[]);
Iy = reshape(Iy,1,1,[]);
Iz = reshape(Iz,1,1,[]);
A = reshape(A,1,1,[]);
L = reshape(L,1,1,[]);
% ================= Cálculos para aumentar a velocidade ================= %
EA_L = E.*A./L;
EIz_L = E.*Iz./L;
EIz_L2 = EIz_L./L;
EIz_L3 =EIz_L2./L;
EIy_L = E.*Iy./L;
EIy_L2 = EIy_L./L;
EIy_L3 =EIy_L2./L;
GIx_L = G.*Ix./L;
Z0 = zeros(1,1,size(EA_L,3));
% ================ Definição da matriz de rigidez local ================= %
K = [ EA_L  Z0         Z0         Z0     Z0        Z0        -EA_L Z0         Z0         Z0     Z0        Z0 ;
      Z0    12*EIz_L3  Z0         Z0     Z0        6*EIz_L2  Z0    -12*EIz_L3 Z0         Z0     Z0        6*EIz_L2;
      Z0    Z0         12*EIy_L3  Z0     -6*EIy_L2 Z0        Z0    Z0         -12*EIy_L3 Z0     -6*EIy_L2 Z0;
      Z0    Z0         Z0         GIx_L  Z0        Z0        Z0    Z0         Z0         -GIx_L Z0        Z0;
      Z0    Z0         -6*EIy_L2  Z0     4*EIy_L   Z0        Z0    Z0         6*EIy_L2   Z0     2*EIy_L   Z0;
      Z0    6*EIz_L2   Z0         Z0     Z0        4*EIz_L   Z0    -6*EIz_L2  Z0         Z0     Z0        2*EIz_L;
      -EA_L Z0         Z0         Z0     Z0        Z0        EA_L  Z0         Z0         Z0     Z0        Z0;
      Z0    -12*EIz_L3 Z0         Z0     Z0        -6*EIz_L2 Z0    12*EIz_L3  Z0         Z0     Z0        -6*EIz_L2;
      Z0    Z0         -12*EIy_L3 Z0     6*EIy_L2  Z0        Z0    Z0         12*EIy_L3  Z0     6*EIy_L2  Z0;
      Z0    Z0         Z0         -GIx_L Z0        Z0        Z0    Z0         Z0         GIx_L  Z0        Z0;
      Z0    Z0         -6*EIy_L2  Z0     2*EIy_L   Z0        Z0    Z0         6*EIy_L2   Z0     4*EIy_L   Z0;
      Z0    6*EIz_L2   Z0         Z0     Z0        2*EIz_L   Z0    -6*EIz_L2  Z0         Z0     Z0        4*EIz_L];

end

% K = [ EA_L  Z0         Z0         Z0     Z0        Z0        -EA_L Z0         Z0         Z0     Z0        Z0 ;
%       Z0    12*EIz_L3  Z0         Z0     Z0        6*EIz_L2  Z0    -12*EIz_L3 Z0         Z0     Z0        6*EIz_L2;
%       Z0    Z0         12*EIy_L3  Z0     -6*EIy_L2 Z0        Z0    Z0         -12*EIy_L3 Z0     -6*EIy_L2 Z0;
%       Z0    Z0         Z0         GIx_L  Z0        Z0        Z0    Z0         Z0         -GIx_L Z0        Z0;
%       Z0    Z0         -6*EIy_L2  Z0     4*EIy_L   Z0        Z0    Z0         6*EIy_L2   Z0     2*EIy_L   Z0;
%       Z0    6*EIz_L2   Z0         Z0     Z0        4*EIz_L   Z0    -6*EIz_L2  Z0         Z0     Z0        2*EIz_L;
%       -EA_L Z0         Z0         Z0     Z0        Z0        EA_L  Z0         Z0         Z0     Z0        Z0;
%       Z0    -12*EIz_L3 Z0         Z0     Z0        -6*EIz_L2 Z0    12*EIz_L3  Z0         Z0     Z0        -6*EIz_L2;
%       Z0    Z0         -12*EIy_L3 Z0     6*EIy_L2  Z0        Z0    Z0         12*EIy_L3  Z0     6*EIy_L2  Z0;
%       Z0    Z0         Z0         -GIx_L Z0        Z0        Z0    Z0         Z0         GIx_L  Z0        Z0;
%       Z0    Z0         -6*EIy_L2  Z0     2*EIy_L   Z0        Z0    Z0         6*EIy_L2   Z0     4*EIy_L   Z0;
%       Z0    6*EIz_L2   Z0         Z0     Z0        2*EIz_L   Z0    -6*EIz_L2  Z0         Z0     Z0        4*EIz_L];

