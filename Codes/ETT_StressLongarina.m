%% Função ETT_StressLongarina() %%
function [ L ] = ETT_StressLongarina( L )

for i = 1:length(L)
    
    for sec = 1:length(L(i).config)
        idx = L(i).posLong(sec):L(i).posLong(sec+1)-1;
        pontos = L(i).pontos(idx,:);
        A = L(i).A(idx); 
        Ix = L(i).Ix(idx); 
        Iy = L(i).Iy(idx); 
        Iz = L(i).Iz(idx);
        Iyz = L(i).Iyz(idx);
        Ex = L(i).Ex(idx); 
        Ey = L(i).Ey(idx); 
        Ez = L(i).Ez(idx); 
        Gx = L(i).Gx(idx); 
        Gy = L(i).Gy(idx); 
        Gz = L(i).Gz(idx);
        
        F = L(i).F(idx,:,:);
        M = L(i).M(idx,:,:);
        
        % ========================== Sanduíche ========================== %
        if strcmp(L(i).config{sec},'sanduiche')
            A_c = L(i).A_c(idx); b = L(i).bc(idx); c = L(i).c(idx); t = L(i).t(idx); h = L(i).h(idx); gama = L(i).gama(idx); Gc = L(i).Gc(idx);
            
            [ Sx_f , T_f , Sx_c , T_c ] = StressSanduiche( A , Ix , Iy , Iz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , A_c , b , c , t , h , gama , Gc );
            
            L(i).Sx(:,idx,:) = Sx_f; 
            L(i).T(:,idx,:) = T_f; 
            L(i).Sx_c(:,idx,:) = Sx_c; 
            L(i).T_c(:,idx,:) = T_c;
            
        elseif strcmp(L(i).config{sec},'sanduiche_I') || strcmp(L(i).config{sec},'sanduiche_H')
            Dvy = L(i).Dvy(idx); Dvz = L(i).Dvz(idx); GJ = L(i).GJ(idx);
            Bv = L(i).Bv(idx); Hv = L(i).Hv(idx); Bn = L(i).Bn(idx); Hn = L(i).Hn(idx);
            En = L(i).En(idx); Ef = L(i).Ef(idx); Gn = L(i).Gn(idx); Gf = L(i).Gf(idx);
            [ Sx_f , T_f , Sx_c , T_c ] = StressSanduiche2( A , Ix , Iy , Iz , Ex , Ey , Ez , En , Ef , Gx , Gy , Gz , Gn , Gf , F , M , Dvy , Dvz , GJ , Bv , Hv , Bn , Hn );
            
            L(i).Sx(:,idx,:) = Sx_f; 
            L(i).T(:,idx,:) = T_f; 
            L(i).Sx_c(:,idx,:) = Sx_c; 
            L(i).T_c(:,idx,:) = T_c;
            
        elseif strcmp(L(i).config{sec},'circular') || strcmp(L(i).config{sec},'tubinho') || strcmp(L(i).config{sec},'cone')
            D_int = L(i).D_int(idx); D_ext = L(i).D_ext(idx);
            
            [ Sx , T ] = StressCircular( F , M , D_int , D_ext , A , Ix , Iy , Iz );
            
            L(i).Sx(:,idx,:) = Sx; 
            L(i).T(:,idx,:) = T; 
            
        elseif strcmp(L(i).config{sec},'retangular')
            b = L(i).bc(idx); c = L(i).c(idx);
            
            [ Sx , T ] = StressRet( A , Ix , Iy , Iz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , b , c );
            
            L(i).Sx(:,idx,:) = Sx; 
            L(i).T(:,idx,:) = T; 
         
        elseif strcmp(L(i).config{sec},'BA')
            Malha_Long = L(i).Malha_Long;
            X = Malha_Long.X(:,idx);
            Z_sup = Malha_Long.Z_sup(:,idx);
            Z_inf = Malha_Long.Z_inf(:,idx);
            t = L(i).t(idx);
            [ Sx , T ] = StressBA( t, A , Ix , Iy , Iz , Iyz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , X , Z_sup , Z_inf , pontos );
            
            L(i).Sx(:,idx,:) = Sx; 
            L(i).T(:,idx,:) = T;
            
        end
    end    
end

end

%% Função ETT_StressBA() %%
function [Sy, T] = StressBA( t, A , Ix , Iy , Iz , Iyz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , X , Z_sup , Z_inf , pontos )
num_pontos = length(X);
num_condicoes = size(F,3);
num_secao = 6;
% ======================== OPERAÇÕES MATRICIAIS ========================= %
P = repmat(permute(F(:,1,:),[2 1 3]),num_secao,1,1);
Vy = repmat(permute(F(:,2,:),[2 1 3]),num_secao,1,1);
Vz = repmat(permute(F(:,3,:),[2 1 3]),num_secao,1,1);

My = repmat(permute(M(:,2,:),[2 1 3]),num_secao,1,1);
Mx = repmat(permute(M(:,1,:),[2 1 3]),num_secao,1,1);
Mz = repmat(permute(M(:,3,:),[2 1 3]),num_secao,1,1);

t_m = repmat(t,num_secao,1,num_condicoes);
A_m = repmat(A,num_secao,1,num_condicoes);
Ix_m = repmat(Ix,num_secao,1,num_condicoes);
Iy_m = repmat(Iy,num_secao,1,num_condicoes);
Iz_m = repmat(Iz,num_secao,1,num_condicoes);
Iyz_m = repmat(Iyz,num_secao,1,num_condicoes);

[m,n] = size(X);
% ============================= COORDENADAS  ============================ %
aux = 1;
for i = 2:num_secao/2
    aux(i) = aux(i-1) + round(m/(num_secao/2));
end

X_secao = X(aux,:) - pontos(:,1)';
Z_sup_secao = flip(Z_sup(aux,:) - pontos(:,3)');
Z_inf_secao = Z_inf(aux,:) - pontos(:,3)';
y_secao = [flip(X_secao);X_secao];
z_secao = [Z_sup_secao;Z_inf_secao];

y = repmat(y_secao,1,1,num_condicoes);
z = repmat(z_secao,1,1,num_condicoes);

% ===================== Primeiro Momento de Inércia ===================== %
dQy = [];
dQz = [];

for i = 1:n
    qy = [];
    qz = [];
    Xq = [flip(X(:,i));X(2:end,i);X(end,i)]';
    Z = [flip(Z_sup(:,i));Z_inf(2:end,i);Z_inf(end,i)]';
    for j = 1:2*m-1
        dA = t(i)*sqrt((Xq(j+1)-Xq(j))^2 + (Z(j+1)-Z(j))^2);
        qz = [qz dA*((Z(j+1)+Z(j))/2 - pontos(i,3))];
        qy = [qy dA*((Xq(j+1)+Xq(j))/2 - pontos(i,1))];
    end
    dQy = [dQy;qy];
    dQz = [dQz;qz];
end
Qy = cumsum(dQy)';
Qz = cumsum(dQz)';


Qy_secao = [Qy(aux,:);Qy((aux+m-3),:)];
Qz_secao = [Qz(aux,:);Qz((aux+m-3),:)];

Qy = repmat(Qy_secao,1,1,num_condicoes);
Qz = repmat(Qz_secao,1,1,num_condicoes);
% =============================== SIGMA Y =============================== %
% Sy = +P./A_m - (My./Iy_m).*z + (Mz./Iz_m).*y;
Sy = +P./A_m - My.*(Iy_m.*z - Iyz_m.*y)./(Iy_m.*Iz_m - Iyz_m.^2) + Mz.*(Iz_m.*y - Iyz_m.*z)./(Iy_m.*Iz_m - Iyz_m.^2);
% ================== Fluxo de Cisalhamento - Cortante =================== %
factor_inercialY = (Vy.*Iy_m - Vz.*Iyz_m)./(Iy.*Iz-Iyz.^2);
factor_inercialZ = (Vz.*Iz_m - Vy.*Iyz_m)./(Iy.*Iz-Iyz.^2);

q_cis = -factor_inercialY.*Qy - factor_inercialZ.*Qz;

T = q_cis./t_m + Mx./(A.*t_m);


end

%% Função ETT_StressCircular() %%
function [Sy,T] = StressCircular(F,M,D_int,D_ext,A,Ix,Iy,Iz)
num_pontos = length(D_int);
num_condicoes = size(F,3);
num_angulos = 6;
% ======================== OPERAÇÕES MATRICIAIS ========================= %
P = repmat(permute(F(:,1,:),[2 1 3]),num_angulos,1,1);
Vy = repmat(permute(F(:,2,:),[2 1 3]),num_angulos,1,1);
Vz = repmat(permute(F(:,3,:),[2 1 3]),num_angulos,1,1);

My = repmat(permute(M(:,2,:),[2 1 3]),num_angulos,1,1);
Mx = repmat(permute(M(:,1,:),[2 1 3]),num_angulos,1,1);
Mz = repmat(permute(M(:,3,:),[2 1 3]),num_angulos,1,1);

A_m = repmat(A,num_angulos,1,num_condicoes);
Ix_m = repmat(Ix,num_angulos,1,num_condicoes);
Iy_m = repmat(Iy,num_angulos,1,num_condicoes);
Iz_m = repmat(Iz,num_angulos,1,num_condicoes);

r1 = repmat(D_int/2,num_angulos,1,num_condicoes);
r2 = repmat(D_ext/2,num_angulos,1,num_condicoes);
% ======================== CÁLCULOS PRELIMINARES ======================== %
razao = (r1.^2 + r1.*r2 + r2.^2)./(r1.^2 + r2.^2);
% ========================= COORDENADAS POLARES ========================= %
theta = repmat(linspace(0,2*pi,num_angulos)',1,num_pontos,num_condicoes);
defasagem = pi/2;
y = cos(theta+defasagem) .* r2;
z = sin(theta+defasagem) .* r2;
% =============================== SIGMA Y =============================== %
Sy = +P./A_m - (My./Iy_m).*z + (Mz./Iz_m).*y;
% =============================== TAU YX ================================ %
T_torsor = (Mx./Ix_m).*r2;
Tyx = T_torsor.*cos(theta) + (4/3)*(Vy./A_m).*razao.*(1 - y.^2./r2.^2);
% =============================== TAU YZ ================================ %
Tyz = T_torsor.*sin(theta) + (4/3)*(Vz./A_m).*razao.*(1 - z.^2./r2.^2);

T = sqrt(Tyx.^2 + Tyz.^2);
end

%% Função StressSanduiche() %%
function [ Sx_f , T_f , Sx_c , T_c ] = StressSanduiche( A , Ix , Iy , Iz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , A_c , b , c , t , h , gama , Gc )

num_elementos = size(A,2);
num_condicoes = size(F,3);
num_secao = 6;

% ======================== Operações Matriciais ========================= %
A_m = repmat(A,num_secao,1,num_condicoes);

Ix_m = repmat(Ix,num_secao,1,num_condicoes);
Iy_m = repmat(Iy,num_secao,1,num_condicoes);
Iz_m = repmat(Iz,num_secao,1,num_condicoes);

A_c_m = repmat(A_c,num_secao,1,num_condicoes);
b_m = repmat(b,num_secao,1,num_condicoes);
h_m = repmat(h,num_secao,1,num_condicoes);
c_m = repmat(c,num_secao,1,num_condicoes);
gama_m = repmat(gama,num_secao,1,num_condicoes);
Gc_m = repmat(Gc,num_secao,1,num_condicoes);

P = repmat(permute(F(:,1,:),[2 1 3]),num_secao,1,1);
Vy = repmat(permute(F(:,2,:),[2 1 3]),num_secao,1,1);
Vz = repmat(permute(F(:,3,:),[2 1 3]),num_secao,1,1);

Mx = repmat(permute(M(:,1,:),[2 1 3]),num_secao,1,1);
My = repmat(permute(M(:,2,:),[2 1 3]),num_secao,1,1);
Mz = repmat(permute(M(:,3,:),[2 1 3]),num_secao,1,1);

% ============================= Coordenadas ============================= %
y = repmat([1;1;0;-1;-1;0],1,num_elementos,num_condicoes) .* b_m/2;
z_f = repmat([-1;1;1;1;-1;-1],1,num_elementos,num_condicoes) .* h_m/2;
% z_c = repmat([-1;1;1;1;-1;-1],1,num_pontos,num_condicoes) .* c_m/2;

% =============================== Sigma X =============================== %
Sx_f = +P./A_m + (My./Iy_m).*z_f - (Mz./Iz_m).*y;                           % as inércias consideradas são só as das faces
Sx_c = zeros(size(Sx_f));                                                   % o core não resiste a tração (Ef >> Ec)

% ================================= Tau ================================= %
K = Ix_m .* Gx;
theta_L = Mx./K;
gb = gama_m.*b;
gy = gama_m.*(y + b_m/2);

T_f = theta_L .* Gx .* (h + t) .* ( sinh(gy).*(cosh(gb)-1)./sinh(gb) - cosh(gy) + 1 );
T_c = theta_L .* Gx .* (h + t) .* sqrt(2.*t.*Gc_m./(h.*Gx)) .* ( cosh(gy).*(cosh(gb)-1)./sinh(gb) - sinh(gy) );

T_c = T_c + (Vy + Vz)./A_c_m;                                                % só o core suporta as tensões de cisalhamento por cortante

end

%% Função StressSanduiche() %%
function [ Sx_f , T_f , Sx_c , T_c ] = StressSanduiche2( A , Ix , Iy , Iz , Ex , Ey , Ez , En , Ef , Gx , Gy , Gz , Gn , Gf , F , M , Dvy , Dvz , GJ , Bv , Hv , Bn , Hn )

num_elementos = size(A,2);
num_condicoes = size(F,3);
num_secao = 8;

% ======================== Operações Matriciais ========================= %
A_m = repmat(A,num_secao,1,num_condicoes);

En_m = repmat(En,num_secao,1,num_condicoes);
Ef_m = repmat(Ef,num_secao,1,num_condicoes);
Gn_m = repmat(Gn,num_secao,1,num_condicoes);
Gf_m = repmat(Gf,num_secao,1,num_condicoes);

Ix_m = repmat(Ix,num_secao,1,num_condicoes);
Iy_m = repmat(Iy,num_secao,1,num_condicoes);
Iz_m = repmat(Iz,num_secao,1,num_condicoes);

Dvy_m = repmat(Dvy,num_secao,1,num_condicoes);
Dvz_m = repmat(Dvz,num_secao,1,num_condicoes);
GJ_m = repmat(GJ,num_secao,1,num_condicoes);

Bv_m = repmat(Bv,num_secao,1,num_condicoes);
Hv_m = repmat(Hv,num_secao,1,num_condicoes);

Bn_m = repmat(Bn,num_secao,1,num_condicoes);
Hn_m = repmat(Hn,num_secao,1,num_condicoes);

P = repmat(permute(F(:,1,:),[2 1 3]),num_secao,1,1);
Vy = repmat(permute(F(:,2,:),[2 1 3]),num_secao,1,1);
Vz = repmat(permute(F(:,3,:),[2 1 3]),num_secao,1,1);

Mx = repmat(permute(M(:,1,:),[2 1 3]),num_secao,1,1);
My = repmat(permute(M(:,2,:),[2 1 3]),num_secao,1,1);
Mz = repmat(permute(M(:,3,:),[2 1 3]),num_secao,1,1);

% ============================= Coordenadas ============================= %
y_n = repmat([1;1;0;-1;-1;-1;0;1],1,num_elementos,num_condicoes) .* Bn_m/2;
y_f = repmat([1;1;0;-1;-1;-1;0;1],1,num_elementos,num_condicoes) .* Bv_m/2;
z_n = repmat([0;1;1;1;0;-1;-1;-1],1,num_elementos,num_condicoes) .* Hn_m/2;
z_f = repmat([0;1;1;1;0;-1;-1;-1],1,num_elementos,num_condicoes) .* Hv_m/2;
r_n = sqrt(y_n.^2 + z_n.^2);
r_f = sqrt(y_f.^2 + z_f.^2);
theta_n = atan2(z_n,y_n);
theta_f = atan2(z_f,y_f);

% =============================== Sigma X =============================== %
Sx_f = +P./A_m + (My.*Ef_m./Dvy_m).*z_f - (Mz.*Ef_m./Dvz_m).*y_f;             
Sx_c = +P./A_m + (My.*En_m./Dvy_m).*z_n - (Mz.*En_m./Dvz_m).*y_n;             

% ===================== Primeiro Momento de Area ======================== %
Qy_nn = Bn_m.*(Hn_m/2 - z_n).*(Hn_m/2 + z_n)/2;
Qy_nf = Bv_m.*(Hv_m/2 - z_n).*(Hv_m/2 + z_n)/2 - Qy_nn;
Qz_nn = Hn_m.*(Bn_m/2 - y_n).*(Bn_m/2 + y_n)/2;
Qz_nf = Hv_m.*(Bv_m/2 - y_n).*(Bv_m/2 + y_n)/2 - Qz_nn;

Qy_fn = Bn_m.*max((Hn_m/2 - z_f),0).*max((Hn_m/2 + z_f),0)/2;
Qy_ff = Bv_m.*(Hv_m/2 - z_f).*(Hv_m/2 + z_f)/2 - Qy_fn;
Qz_fn = Hn_m.*max((Bn_m/2 - y_f),0).*max((Bn_m/2 + y_f),0)/2;
Qz_ff = Hv_m.*(Bv_m/2 - y_f).*(Bv_m/2 + y_f)/2 - Qz_fn;

% =============== Primeiro Momento de Area Equivalente ================== %
ESz_n = En_m.*Qy_nn + Ef_m.*Qy_nf;
ESy_n = En_m.*Qz_nn + Ef_m.*Qz_nf;
ESz_f = En_m.*Qy_fn + Ef_m.*Qy_ff;
ESy_f = En_m.*Qz_fn + Ef_m.*Qz_ff;

% ============================= Cortante ================================ %
Cz_c = - Vz.*ESz_n./(Bn_m.*Dvy);
Cy_c = Vy.*ESy_n./(Hn_m.*Dvz);
Cz_f = -Vz.*ESz_f./(Bv_m.*Dvy);
Cy_f = Vy.*ESy_f./(Hv_m.*Dvz);

% ============================== Torsor ================================= %
Tsz_c = -cos(theta_n).*Mx.*Gn_m.*r_n./GJ_m;
Tsy_c = -sin(theta_n).*Mx.*Gn_m.*r_n./GJ_m;
Tsz_f = -cos(theta_f).*Mx.*Gf_m.*r_f./GJ_m;
Tsy_f = -sin(theta_f).*Mx.*Gf_m.*r_f./GJ_m;

% ================================ Tau ================================== %
T_c = sqrt((Cz_c + Tsz_c).^2 + (Cy_c + Tsy_c).^2);
T_f = sqrt((Cz_f + Tsz_f).^2 + (Cy_f + Tsy_f).^2);
% K = Ix_m .* Gx;
% theta_L = Mx./K;
% gb = gama_m.*b;
% gy = gama_m.*(y + b_m/2);
% 
% T_f = theta_L .* Gx .* (h + t) .* ( sinh(gy).*(cosh(gb)-1)./sinh(gb) - cosh(gy) + 1 );
% T_c = theta_L .* Gx .* (h + t) .* sqrt(2.*t.*Gc_m./(h.*Gx)) .* ( cosh(gy).*(cosh(gb)-1)./sinh(gb) - sinh(gy) );
% 
% T_c = T_c + (Vy + Vz)./A_c_m;                                                % só o core suporta as tensões de cisalhamento por cortante

end

%% Função StressSanduiche() %%

function [ Sx , T ] = StressRet( A , Ix , Iy , Iz , Ex , Ey , Ez , Gx , Gy , Gz , F , M , b , c )

num_elementos = size(A,2);
num_condicoes = size(F,3);
num_secao = 6;

% ======================== Operações Matriciais ========================= %
A_m = repmat(A,num_secao,1,num_condicoes);

Ix_m = repmat(Ix,num_secao,1,num_condicoes);
Iy_m = repmat(Iy,num_secao,1,num_condicoes);
Iz_m = repmat(Iz,num_secao,1,num_condicoes);

b_m = repmat(b,num_secao,1,num_condicoes);
c_m = repmat(c,num_secao,1,num_condicoes);

P = repmat(permute(F(:,1,:),[2 1 3]),num_secao,1,1);
Vy = repmat(permute(F(:,2,:),[2 1 3]),num_secao,1,1);
Vz = repmat(permute(F(:,3,:),[2 1 3]),num_secao,1,1);

Mx = repmat(permute(M(:,1,:),[2 1 3]),num_secao,1,1);
My = repmat(permute(M(:,2,:),[2 1 3]),num_secao,1,1);
Mz = repmat(permute(M(:,3,:),[2 1 3]),num_secao,1,1);

% ============================= Coordenadas ============================= %
y = repmat([1;1;0;-1;-1;0],1,num_elementos,num_condicoes) .* b_m/2;
z_f = repmat([-1;1;1;1;-1;-1],1,num_elementos,num_condicoes) .* c_m/2;
% z_c = repmat([-1;1;1;1;-1;-1],1,num_pontos,num_condicoes) .* c_m/2;

% =============================== Sigma X =============================== %
Sx_f = +P./A_m + (My./Iy_m).*z_f - (Mz./Iz_m).*y;                           % as inércias consideradas são só as das faces
Sx_c = zeros(size(Sx_f));                                                   % o core não resiste a tração (Ef >> Ec)

% ================================= Tau ================================= %
K = Ix_m .* Gx;
theta_L = Mx./K;
gb = gama_m.*b;
gy = gama_m.*(y + b_m/2);

T_f = theta_L .* Gx .* (h + t) .* ( sinh(gy).*(cosh(gb)-1)./sinh(gb) - cosh(gy) + 1 );
T_c = theta_L .* Gx .* (h + t) .* sqrt(2.*t.*Gc_m./(h.*Gx)) .* ( cosh(gy).*(cosh(gb)-1)./sinh(gb) - sinh(gy) );

T_c = T_c + (Vy + Vz)./A_c_m;                                                % só o core suporta as tensões de cisalhamento por cortante

end