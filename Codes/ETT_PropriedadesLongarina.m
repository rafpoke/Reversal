function [Longarina] = ETT_PropriedadesLongarina( Longarina )
%% CARREGANDO VALORES %%
config = {Longarina.config};

%% CHAMANDO FUNÇÕES ESPECÍFICAS %%

Longarina = inicializaLongarina( Longarina );

for i = 1:length(Longarina)
    for j = 1:length(Longarina(i).config)
        if strcmp(Longarina(i).config{j},'sanduiche')
            Longarina = propSand( Longarina , i , j );
        elseif strcmp(Longarina(i).config{j},'tubinho')
            Longarina = propTubinho( Longarina , i , j );
        elseif strcmp(Longarina(i).config{j},'circular')
            Longarina = propCircular( Longarina , i , j );
        elseif strcmp(Longarina(i).config{j},'cone')
            Longarina = propCone( Longarina, i , j );
        elseif strcmp(Longarina(i).config{j},'retangular')
            Longarina = propRet( Longarina, i , j );
        elseif strcmp(Longarina(i).config{j},'BA')
            Longarina = propBA(Longarina, i , j );
        elseif strcmp(Longarina(i).config{j},'sanduiche_I')
            Longarina = propSandI(Longarina, i , j );
        elseif strcmp(Longarina(i).config{j},'sanduiche_H')
            Longarina = propSandH(Longarina, i , j );            
            
        end
    end
end

end

%% Função propSandI() %%
function [ Longarina ] = propSandI( Longarina , idx_long , idx_sec)
%% Carregando Valores %%
L_w = Longarina(idx_long);
 
H = L_w.H([idx_sec,idx_sec+1]);
Hmax_sec = L_w.Hmax_sec([idx_sec,idx_sec+1]);
n_lam = L_w.n_lam(idx_sec,:);
f_lam = L_w.f_lam(idx_sec);
B = L_w.B([idx_sec,idx_sec+1]);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
material_c = L_w.material_c{idx_sec};
pontos = L_w.pontos(posLong(1):posLong(2),:);

H = Hmax_sec.*H;

Hv = linspaceDsec(H,posLong);
Bv = linspaceDsec(B,posLong);
%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);
 
% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);
Ef = material.E1 * ones(1,num_elementos);
En = material_c.E1 * ones(1,num_elementos);
 
Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);
Gn = material_c.G1 * ones(1,num_elementos);
Gf = material.G1 * ones(1,num_elementos);
 
% e_volta = material.e;
 
% ============================ INICIALIZAÇÃO ============================ %
vet_lam = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));
 
% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_lam;
 
% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_lam(:) = n_lam(2);
vet_lam(antes_troca) = n_lam(1);
 
% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
h = vet_lam * e_volta;
% ========================= VETOR DE DIÂMETROS ========================== %
Hn = Hv - 2*h;
Hm = Hn + h;
 
% =========================== VETOR DE ÁREAS ============================ %
Af = 2*h.*Bv;
An = Hn.*Bv;
A = Af + An;

% ========================== VETOR DE INÉRCIA =========================== %
Iyf = 2*(Bv.*(h.^3)/12 + Bv.*h.*(Hm.^2)/4);
Iyn = Bv.*(Hn.^3)/12;
Iy = Iyf + Iyn;

Izf = 2*h.*(Bv.^3)/12;
Izn = Hn.*(Bv.^3)/12;
Iz = Izf + Izn;

Ixf = Iyf + Iyn;
Ixn = Iyn + Izn;
Ix = Iy + Iz;

Dvy = Ef.*Iyf + En.*Iyn;
Dvz = Ef.*Izf + En.*Izn;

GJ = Gn.*Ixn + Gf.*Ixf;

Iyz = zeros(size(Ix));

% ================================ Massa ================================ %
rhof = material.rho;
rhon = material_c.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* (Af * rhof + An*rhon);
 
%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;
Longarina(idx_long).Ef(idx_elementos) = Ef;
Longarina(idx_long).En(idx_elementos) = En;
 
Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;
Longarina(idx_long).Gn(idx_elementos) = Gn;
Longarina(idx_long).Gf(idx_elementos) = Gf;
 
Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).Dvy(idx_elementos) = Dvy;
Longarina(idx_long).Dvz(idx_elementos) = Dvz;
Longarina(idx_long).GJ(idx_elementos) = GJ;

Longarina(idx_long).A(idx_elementos) = A;
 
Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);

Longarina(idx_long).Hn(idx_elementos) = Hn;
Longarina(idx_long).Bn(idx_elementos) = Bv;
Longarina(idx_long).Hv(idx_elementos) = Hv;
Longarina(idx_long).Bv(idx_elementos) = Bv;
 
%% Função linspaceDsec() %%
% Cria um vetor com os diâmetros entre as seções
function [Dv] = linspaceDsec (D_sec, posLong)
Dv(1) = D_sec(1);
num_posLong = length(posLong);
for i = 1:num_posLong-1
    D1 = linspace(D_sec(i),D_sec(i+1),posLong(i+1)-posLong(i));
    D = D1(2:(posLong(i+1)-posLong(i)));
    Dv = [Dv,D];
end
end
end

%% Função propSandH() %%
function [ Longarina ] = propSandH( Longarina , idx_long , idx_sec)
%% Carregando Valores %%
L_w = Longarina(idx_long);
 
B = L_w.B([idx_sec,idx_sec+1]);
n_lam = L_w.n_lam(idx_sec,:);
f_lam = L_w.f_lam(idx_sec);
H = L_w.H([idx_sec,idx_sec+1]);
Hmax_sec = L_w.Hmax_sec([idx_sec,idx_sec+1]);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
material_c = L_w.material_c{idx_sec};
pontos = L_w.pontos(posLong(1):posLong(2),:);

H = Hmax_sec.*H;

Bv = linspaceDsec(B,posLong);
Hv = linspaceDsec(H,posLong);
 
%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);
 
% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);
Ef = material.E1 * ones(1,num_elementos);
En = material_c.E1 * ones(1,num_elementos);
 
Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);
Gf = material.G1 * ones(1,num_elementos);
Gn = material_c.G1 * ones(1,num_elementos);
 
% e_volta = material.e;
 
% ============================ INICIALIZAÇÃO ============================ %
vet_lam = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));
 
% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_lam;
 
% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_lam(:) = n_lam(2);
vet_lam(antes_troca) = n_lam(1);
 
% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
b = vet_lam * e_volta;
% ========================= VETOR DE DIÂMETROS ========================== %
Bn = Bv - 2*b;
Bm = Bn + b;
 
% =========================== VETOR DE ÁREAS ============================ %
Af = 2*b.*Hv;
An = Bn.*Hv;
A = Af + An;

% ========================== VETOR DE INÉRCIA =========================== %
Iyf = 2*b.*(Hv.^3)/12;
Iyn = Bn.*(Hv.^3)/12;
Iy = Iyf + Iyn;

Izf = 2*(Hv.*(b.^3)/12 + b.*Hv.*(Bm.^2)/4);
Izn = Hv.*(Bn.^3)/12;
Iz = Izf + Izn;

Ixf = Iyf + Iyn;
Ixn = Iyn + Izn;
Ix = Iy + Iz;

Dvy = Ef.*Iyf + En.*Iyn;
Dvz = Ef.*Izf + En.*Izn;

GJ = Gn.*Ixn + Gf.*Ixf;

Iyz = zeros(size(Ix));

% ================================ Massa ================================ %
rhof = material.rho;
rhon = material_c.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* (Af * rhof + An*rhon);
 
%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;
Longarina(idx_long).Ef(idx_elementos) = Ef;
Longarina(idx_long).En(idx_elementos) = En;
 
Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;
Longarina(idx_long).Gf(idx_elementos) = Gf;
Longarina(idx_long).Gn(idx_elementos) = Gn;
 
Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).Dvy(idx_elementos) = Dvy;
Longarina(idx_long).Dvz(idx_elementos) = Dvz;
Longarina(idx_long).GJ(idx_elementos) = GJ;
 
Longarina(idx_long).A(idx_elementos) = A;
 
Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);

Longarina(idx_long).Hn(idx_elementos) = Hv;
Longarina(idx_long).Bn(idx_elementos) = Bn;
Longarina(idx_long).Bv(idx_elementos) = Bv;
Longarina(idx_long).Hv(idx_elementos) = Hv;
 
%% Função linspaceDsec() %%
% Cria um vetor com os diâmetros entre as seções
function [Dv] = linspaceDsec (D_sec, posLong)
Dv(1) = D_sec(1);
num_posLong = length(posLong);
for i = 1:num_posLong-1
    D1 = linspace(D_sec(i),D_sec(i+1),posLong(i+1)-posLong(i));
    D = D1(2:(posLong(i+1)-posLong(i)));
    Dv = [Dv,D];
end
end
end

%% Função propBA() %%
function [ Longarina ] = propBA( Longarina , idx_long , idx_sec )
%% Carregando Valores %%
L_w = Longarina(idx_long);
 
n_voltas = L_w.n_voltas(idx_sec,:);
f_voltas = L_w.f_voltas(idx_sec);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
pontos = L_w.pontos(posLong(1):posLong(2),:);

Malha = struct;
Malha.X = L_w.Malha_Long.X(:,posLong(1):posLong(2));
Malha.Y = L_w.Malha_Long.X(:,posLong(1):posLong(2));
Malha.Z_sup = L_w.Malha_Long.Z_sup(:,posLong(1):posLong(2));
Malha.Z_inf = L_w.Malha_Long.Z_inf(:,posLong(1):posLong(2));


[m,n] = size(Malha.X);
 
%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);
 
% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);
 
Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);
 
% e_volta = material.e;
 
% ============================ INICIALIZAÇÃO ============================ %
vet_voltas = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));
 
% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_voltas;
 
% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_voltas(:) = n_voltas(2);
vet_voltas(antes_troca) = n_voltas(1);
 
% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
t = vet_voltas * e_volta;
 
% =========================== VETOR DE ÁREAS ============================ %
P = zeros(n-1,1);

for i = 1:n-1
    p = 0;
    X = [flip(Malha.X(:,i));Malha.X(2:end,i);Malha.X(end,i)]';
    Z = [flip(Malha.Z_sup(:,i));Malha.Z_inf(2:end,i);Malha.Z_sup(end,i)]';
    for j = 1:2*m-1
        p = sqrt((X(j+1)-X(j))^2 + (Z(j+1)-Z(j))^2) + p;
    end
    P(i) = p;
end

A = P'.*t;
 
% ========================== VETOR DE INÉRCIA =========================== %
Iy = zeros(n-1,1);
Iz = zeros(n-1,1);
Iyz = zeros(n-1,1);

for i = 1:n-1
    iy = 0;
    iz = 0;
    iyz = 0;
    X = [flip(Malha.X(:,i));Malha.X(2:end,i);Malha.X(end,i)]';
    Z = [flip(Malha.Z_sup(:,i));Malha.Z_inf(2:end,i);Malha.Z_sup(end,i)]';
    for j = 1:2*m-1
        iy = t(i)*sqrt((X(j+1)-X(j))^2 + (Z(j+1)-Z(j))^2)*(pontos(i,3)-(Z(j+1)+Z(j))/2)^2 + iy;
        iz = t(i)*sqrt((X(j+1)-X(j))^2 + (Z(j+1)-Z(j))^2)*(pontos(i,1)-(X(j+1)+X(j))/2)^2 + iz;
        iyz = t(i)*sqrt((X(j+1)-X(j))^2 + (Z(j+1)-Z(j))^2)*(pontos(i,1)-(X(j+1)+X(j))/2)*(pontos(i,3)-(Z(j+1)+Z(j))/2) + iyz;
    end
    Iy(i) = iy;
    Iz(i) = iz;
    Iyz(i) = iyz;
end

Ix = Iy + Iz;
% ================================ Massa ================================ %
rho = material.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* A * rho;
 
%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;
 
Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;
 
Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).A(idx_elementos) = A;
Longarina(idx_long).t(idx_elementos) = t;
 
Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);
 
end

%% Função propCircular() %%
function Longarina = propCircular( Longarina , idx_long , idx_sec )
%% Carregando Valores %%
L_w = Longarina(idx_long);

D = L_w.D(idx_sec);
n_voltas = L_w.n_voltas(idx_sec,:);
f_voltas = L_w.f_voltas(idx_sec);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
pontos = L_w.pontos([posLong(1):posLong(2)],:);

%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);

% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);

Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);

% ============================ INICIALIZAÇÃO ============================ %
vet_voltas = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));

% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_voltas;

% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_voltas(:) = n_voltas(2);
vet_voltas(antes_troca) = n_voltas(1);

% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
t = vet_voltas * e_volta;

% ========================= VETOR DE DIÂMETROS ========================== %
D_ext = D + 2*t;
D_int = D * ones(1,num_elementos);

% =========================== VETOR DE ÁREAS ============================ %
A = pi * (D_ext.^2 - D_int.^2)/4;

% ========================== VETOR DE INÉRCIA =========================== %
Iy = (pi/4) * ((D_ext/2).^4 - (D_int/2).^4);
Iz = Iy;
Ix = Iy + Iz;
Iyz = zeros(size(Ix));

% ================================ Massa ================================ %
rho = material.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* A * rho;

%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;

Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;

Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).A(idx_elementos) = A;

Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);

Longarina(idx_long).D_int(idx_elementos) = D_int;
Longarina(idx_long).D_ext(idx_elementos) = D_ext;

end

%% Função propTubinho() %%
function Longarina = propTubinho( Longarina , idx_long , idx_sec )
%% Carregando Valores %%
L_w = Longarina(idx_long);

D_interno = L_w.D_interno(idx_sec);
D_externo = L_w.D_externo(idx_sec);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
pontos = L_w.pontos([posLong(1):posLong(2)],:);

%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);

% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);

Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);

% =============================== Inercias ============================== %
Ix = pi * (D_externo^4 - D_interno^4)/32 * ones(1,num_elementos);
Iy = pi * (D_externo^4 - D_interno^4)/64 * ones(1,num_elementos);
Iz = Iy;
Iyz = zeros(size(Ix));

% ================================ Áreas ================================ %
A = pi * (D_externo^2 - D_interno^2)/4 * ones(1,num_elementos);

% ========================= VETOR DE DIÂMETROS ========================== %
D_ext = D_externo * ones(1,num_elementos);
D_int = D_interno * ones(1,num_elementos);

% ================================ Massa ================================ %
rho = material.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* A * rho;

%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;

Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;

Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).A(idx_elementos) = A;

Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);

Longarina(idx_long).D_int(idx_elementos) = D_int;
Longarina(idx_long).D_ext(idx_elementos) = D_ext;


end

%% Função propSand() %%
function [Longarina] = propSand( Longarina , idx_long , idx_sec )
%% CARREGANDO VALORES %%
L_w = Longarina(idx_long);

D_I = L_w.D_I{idx_sec};
beta = L_w.beta{idx_sec};
b = L_w.b{idx_sec};
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
material_c = L_w.material_c{idx_sec};
pontos = L_w.pontos([posLong(1):posLong(2)],:);

%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);

% ==================== Definindo variáveis contínuas ==================== %
x = linspace(0,1,num_elementos);

D_I = polyval( F_Spline(D_I(1),D_I(2),D_I(3)) , x );
beta = polyval( F_Spline(beta(1),beta(2),beta(3)) , x );
b = polyval( F_Spline(b(1),b(2),b(3)) , x );

% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);

Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);

Gc = material_c.G1 * ones(1,num_elementos);

% ========================== Rigidez torsional ========================== %
K_I = beta .* D_I .* (Gx./Ex);

% =========================== Core thickness ============================ %
rho_f = material.rho;
rho_c = material_c.rho;

c_opt1 = 2 * (rho_f .* K_I ./ (rho_c .* Gx)).^(1/3);
c_opt2 = 2 * (rho_f .* D_I ./ (rho_c .* Ex)).^(1/3);

ratio = K_I./D_I - Gx./Ex;
c_opt = c_opt1 .* (ratio >= 0) + c_opt2 .* (ratio < 0);

% c_opt = round(c_opt * 1000)/1000;                                           % aproxima c_opt para o inteiro mais próximo (em mm)
c_opt(c_opt < 1e-3) = 1e-3;

% ============================ Face thickness =========================== %
t_opt = rho_c * c_opt / (4*rho_f);

t_opt = round(t_opt * 4000)/4000;                                           % aproxima t_opt para o número de camadas mais próximo (passo de 0.25 mm)
t_opt(t_opt < 0.25e-3) = 0.25e-3;

% ======================== Recalculando rigidez ========================= %
D_I = Ex .* t_opt .* c_opt.^2/2;
K_I = Gx .* t_opt .* c_opt.^2/2;

% =============================== Inercias ============================== %
Ix3 = (K_I./Gx) .* b;  % não está concordando com o elementos finitos
Iy = (D_I./Ex) .* b;
Iz = t_opt .* c_opt.^2/2 .* b;  % conferir se isso tá certo (acho que não)
Iyz = zeros(size(Iz));
% ================================ Áreas ================================ %
A_f = 2 * b .* t_opt;
A_c = b .* c_opt;

% ================================ Massa ================================ %
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* (A_f * rho_f + A_c * rho_c);

% ================================== h ================================== %
h = c_opt + 2 * t_opt;

% ================================ gama ================================= %
gama = sqrt( (2 .* h .* Gc)./(t_opt .* Gx .* (h + t_opt).^2) );

gb = gama .* b;
Ix = 2.*t_opt.*(h+t_opt).^2.*b .* ( 1 - 2.*(cosh( gb ) - 1)./( gb .* sinh( gb )));    % fórmula alternativa que parece melhor


%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;

Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;

Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;

Longarina(idx_long).A(idx_elementos) = A_f;

Longarina(idx_long).A_c(idx_elementos) = A_c;
Longarina(idx_long).D_Ic(idx_elementos) = D_I;
Longarina(idx_long).K_Ic(idx_elementos) = K_I;

Longarina(idx_long).t(idx_elementos) = t_opt;
Longarina(idx_long).c(idx_elementos) = c_opt;
Longarina(idx_long).bc(idx_elementos) = b;
Longarina(idx_long).h(idx_elementos) = h;
Longarina(idx_long).gama(idx_elementos) = gama;
Longarina(idx_long).Gc(idx_elementos) = Gc;

Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);

Longarina(idx_long).Ix3(idx_elementos) = Ix3;


end

%% Função propCone() %%
function [ Longarina ] = propCone( Longarina , idx_long , idx_sec )
%% Carregando Valores %%
L_w = Longarina(idx_long);
 
D_sec = L_w.D_sec([idx_sec,idx_sec+1]);
n_voltas = L_w.n_voltas(idx_sec,:);
f_voltas = L_w.f_voltas(idx_sec);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
pontos = L_w.pontos([posLong(1):posLong(2)],:);
 
 
Dv = linspaceDsec(D_sec,posLong);
 
%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);
 
% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);
 
Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);
 
% e_volta = material.e;
 
% ============================ INICIALIZAÇÃO ============================ %
vet_voltas = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));
 
% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_voltas;
 
% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_voltas(:) = n_voltas(2);
vet_voltas(antes_troca) = n_voltas(1);
 
% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
t = vet_voltas * e_volta;
% ========================= VETOR DE DIÂMETROS ========================== %
D_ext = Dv + 2*t;
D_int = Dv .* ones(1,num_elementos);
 
% =========================== VETOR DE ÁREAS ============================ %
A = pi * (D_ext.^2 - Dv.^2)/4;
 
% ========================== VETOR DE INÉRCIA =========================== %
Iy = (pi/4) * ((D_ext/2).^4 - (Dv/2).^4);
Iz = Iy;
Ix = Iy + Iz;
Iyz = zeros(size(Ix));
 
% ================================ Massa ================================ %
rho = material.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* A * rho;
 
%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;
 
Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;
 
Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;
 
Longarina(idx_long).A(idx_elementos) = A;
 
Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);
 
Longarina(idx_long).D_int(idx_elementos) = D_int;
Longarina(idx_long).D_ext(idx_elementos) = D_ext;
 
%% Função linspaceDsec() %%
% Cria um vetor com os diâmetros entre as seções
function [Dv] = linspaceDsec (D_sec, posLong)
Dv(1) = D_sec(1);
num_posLong = length(posLong);
for i = 1:num_posLong-1
    D1 = linspace(D_sec(i),D_sec(i+1),posLong(i+1)-posLong(i));
    D = D1(2:(posLong(i+1)-posLong(i)));
    Dv = [Dv,D];
end
end
end

%% Função propCone() %%
function [ Longarina ] = propRet( Longarina , idx_long , idx_sec )
%% Carregando Valores %%
L_w = Longarina(idx_long);
 
base_sec = L_w.base_sec([idx_sec,idx_sec+1]);
alma_sec = L_w.alma_sec([idx_sec,idx_sec+1]);
n_voltas = L_w.n_voltas(idx_sec,:);
f_voltas = L_w.f_voltas(idx_sec);
posLong = L_w.posLong([idx_sec,idx_sec+1]);
material = L_w.material{idx_sec};
pontos = L_w.pontos([posLong(1):posLong(2)],:);
 
  
%% DEFININDO PROPRIEDADES %%
idx_elementos = posLong(1):(posLong(2)-1);
num_elementos = length(idx_elementos);
 
% ================================ E / G ================================ %
Ex = material.E1 * ones(1,num_elementos);
Ey = material.E2 * ones(1,num_elementos);
Ez = material.E2 * ones(1,num_elementos);
 
Gx = material.G1 * ones(1,num_elementos);
Gy = material.G2 * ones(1,num_elementos);
Gz = material.G2 * ones(1,num_elementos);
 
% e_volta = material.e;
 
% ============================ INICIALIZAÇÃO ============================ %
vet_voltas = zeros(1,num_elementos);
pontos_elem = 1/2 * (pontos(1:end-1,:) + pontos(2:end,:));
 
% ====================== POSIÇÃO DE TROCA DE VOLTAS ===================== %
b_long = pontos(end,2) - pontos(1,2);
y_troca = pontos(1,2) + b_long * f_voltas;
 
% =========================== VETOR DE VOLTAS =========================== %
antes_troca = pontos_elem(:,2) <= y_troca;
vet_voltas(:) = n_voltas(2);
vet_voltas(antes_troca) = n_voltas(1);
 
% ========================= VETOR DE ESPESSURA ========================== %
e_volta = 0.25e-3;
t = vet_voltas * e_volta;
% ========================= VETOR DE DIÂMETROS ========================== %
Bv = linspace(base_sec(1),base_sec(2),num_elementos);
Av = linspace(alma_sec(1),alma_sec(2),num_elementos);

A_ext = Av + 2*t;
A_int = Av;
B_ext = Bv + 2*t;
B_int = Bv;

% =========================== VETOR DE ÁREAS ============================ %
A = A_ext.*B_ext - A_int.*B_int;
 
% ========================== VETOR DE INÉRCIA =========================== %
Iy = 1/12 * (A_ext.^3 + B_ext) - 1/12 * (A_int.^3 + B_int);
Iz = 1/12 * (A_ext + B_ext.^3) - 1/12 * (A_int + B_int.^3);
Ix = Iy + Iz;
Iyz = zeros(size(Ix));
% ================================ Massa ================================ %
rho = material.rho;
dx = sqrt(sum(diff(pontos).^2,2));
m = dx' .* A * rho;
 
%% OUTPUTS %%
Longarina(idx_long).Ex(idx_elementos) = Ex;
Longarina(idx_long).Ey(idx_elementos) = Ey;
Longarina(idx_long).Ez(idx_elementos) = Ez;
 
Longarina(idx_long).Gx(idx_elementos) = Gx;
Longarina(idx_long).Gy(idx_elementos) = Gy;
Longarina(idx_long).Gz(idx_elementos) = Gz;
 
Longarina(idx_long).Ix(idx_elementos) = Ix;
Longarina(idx_long).Iy(idx_elementos) = Iy;
Longarina(idx_long).Iz(idx_elementos) = Iz;
Longarina(idx_long).Iyz(idx_elementos) = Iyz;
 
Longarina(idx_long).A(idx_elementos) = A;
 
Longarina(idx_long).m_vet(idx_elementos) = m;
Longarina(idx_long).m = Longarina(idx_long).m + sum(m);
 
Longarina(idx_long).alma(idx_elementos) = A_ext;
Longarina(idx_long).base(idx_elementos) = B_ext;

end

%% Função inicializaLongarina() %%
function [Longarina] = inicializaLongarina( Longarina )

for i = 1:length(Longarina)

    num_pontos = size(Longarina(i).pontos,1);
    num_elementos = num_pontos - 1;

    Longarina(i).A = zeros(1,num_elementos);

    Longarina(i).Ix = zeros(1,num_elementos);
    Longarina(i).Iy = zeros(1,num_elementos);
    Longarina(i).Iz = zeros(1,num_elementos);

    Longarina(i).Ex = zeros(1,num_elementos);
    Longarina(i).Ey = zeros(1,num_elementos);
    Longarina(i).Ez = zeros(1,num_elementos);

    Longarina(i).Gx = zeros(1,num_elementos);
    Longarina(i).Gy = zeros(1,num_elementos);
    Longarina(i).Gz = zeros(1,num_elementos);

    Longarina(i).D_Ic = zeros(1,num_elementos);
    Longarina(i).K_Ic = zeros(1,num_elementos);

    Longarina(i).t = zeros(1,num_elementos);
    Longarina(i).c = zeros(1,num_elementos);
    Longarina(i).bc = zeros(1,num_elementos);
    
    Longarina(i).D_int = zeros(1,num_elementos);
    Longarina(i).D_ext = zeros(1,num_elementos);
    
    Longarina(i).alma = zeros(1,num_elementos);
    Longarina(i).base = zeros(1,num_elementos);

    Longarina(i).m_vet = zeros(1,num_elementos);
    Longarina(i).m = 0;

end

end
