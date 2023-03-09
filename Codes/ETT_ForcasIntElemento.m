%% Função ETT_ForcasIntElemento() %%
% ============================== Descrição ============================== %
% Essa função recebe os deslocamento lineares de cada elemento e calcula as
% forças internas agindo em cada elemento à partir das equações (6.44) do
% Megson. CUIDADO, TEM TERMOS COM SINAL TROCADO NO MEGSON.

function [F_int,M_int] = ETT_ForcasIntElemento( u_loc , elem )
%% Carregando Valores %%
% ==================== Carregando struct_propriedades =================== %
E = elem.E;
G = elem.G;
Ix = elem.Ix;
Iy = elem.Iy;
Iz = elem.Iz;
A = elem.A;
L = elem.L;

%% Corpo da Função %%
% ===================== Deslocamentos nas longarinas ==================== %
% u_loc = u_loc(:,:,idx_elem);
num_elem = size(u_loc,3);
% ============================ Inicialização ============================ %
n_dim = 6;
R_int = zeros(n_dim,1,num_elem);

% =========================== Forças internas =========================== %
for i = 1:num_elem
%     k = idx_elem(i);
    k = i;
    no1 = (1-1)*n_dim;
    no2 = (2-1)*n_dim;
    % ================================ Fx =============================== %
    R_int(1,:,i) = E(k)*A(k)/L(k) * [-1 1] * u_loc([no1+1,no2+1],:,i);
    % ================================ Mx =============================== %
    R_int(4,:,i) = G(k)*Ix(k)/L(k) * [-1 1] * u_loc([no1+4,no2+4],:,i);
    % ================================ Fy =============================== %
    R_int(2,:,i) = E(k)*Iz(k) * [12/L(k)^3 6/L(k)^2 -12/L(k)^3 6/L(k)^2] * u_loc([no1+2,no1+6,no2+2,no2+6],:,i);
    % ================================ Fz =============================== %
    R_int(3,:,i) = E(k)*Iy(k) * [12/L(k)^3 -6/L(k)^2 -12/L(k)^3 -6/L(k)^2] * u_loc([no1+3,no1+5,no2+3,no2+5],:,i);
    % ================================ Mz =============================== %
    m0 = E(k)*Iz(k) * [-6/L(k)^2 , -4/L(k) , 6/L(k)^2 , -2/L(k) ] * u_loc([no1+2,no1+6,no2+2,no2+6],:,i);
    mL = m0 + L(k) * R_int(2,:,i);
    R_int(6,:,i) = maxabs([m0,mL]);
    % ================================ My =============================== %
    m0 = E(k)*Iy(k) * [-6/L(k)^2 , 4/L(k) , 6/L(k)^2 , 2/L(k) ] * u_loc([no1+3,no1+5,no2+3,no2+5],:,i);
    mL = m0 + L(k) * R_int(3,:,i);
    R_int(5,:,i) = maxabs([m0,mL]);
end

% ====================== Mudando formato da matriz ====================== %
% Essa seção do código deixa as forças e os momentos internos no formato
% usual.

F_int = R_int(1:3,:,:);
M_int = R_int(4:end,:,:);

F_int = reshape(permute(F_int,[2 3 1]),[],3);
M_int = reshape(permute(M_int,[2 3 1]),[],3);


end

%% Função maxabs() %%
function max_abs = maxabs(x)

[max_abs,idx] = max(abs(x));
max_abs = max_abs * sign(x(idx));

end
