%% Função ETT_CossenosDiretores() %%
% ============================== Descrição ============================== %
% Essa função recebe os nós e elementos e calcula as matrizes de cossenos
% diretores de cada elemento. Inicialmente é calculado o versor de cada
% elemento e caso ele tenha componente x negativa (< 0) o elemento é
% invertido (caso contrário o MEF não funciona).

function [M_rot,elem] = ETT_CossenosDiretores( nos , elem )
%% Carregando Valores %%
pontos = nos.pontos;
conectividade = elem.conect;

%% Corpo da Função %%
% ================== Definição dos vetores dos elementos ================ %
vet_elem = pontos(conectividade(:,2),:) - pontos(conectividade(:,1),:);
% =================== Arruma sentido da conectividade =================== %
sentido_errado = vet_elem(:,1) < 0;
vet_elem(sentido_errado,:) = -vet_elem(sentido_errado,:);
conectividade(sentido_errado,:) = fliplr(conectividade(sentido_errado,:));
% ================= Definição dos versores dos elementos ================ %
norm_elem = sqrt(vet_elem(:,1).^2 + vet_elem(:,2).^2 + vet_elem(:,3).^2);
vers_elem = bsxfun(@rdivide,vet_elem,norm_elem);
% =================== Definição da base dos elementos =================== %
X_local = vers_elem;
Y_local = zeros(size(X_local));
Z_local = zeros(size(X_local));
Z = [0 0 1];
% A função null(x) recebe um versor e calcula os versores de uma base
% ortonormal positiva que são adotados como os versores do sistema de
% coordenadas local de cada elemento.
% A função null(x) devolve nan se x == Z, então são programados casos
% especiais.
for i = 1:size(X_local,1)
    x = X_local(i,:);
    if x == Z
        Y_local(i,:) = [-1 0 0];
        Z_local(i,:) = [0 1 0];
    elseif x == -Z
        Y_local(i,:) = [1 0 0];
        Z_local(i,:) = [0 -1 0];
    else
        yz=null(x).';
        Y_local(i,:) = yz(1,:);
        Z_local(i,:) = yz(2,:);
    end
end

% ======================== Ângulos de Euler ZYX ========================= %
% [alpha,beta,gama] = EulerAngles(Y_local,Z_local);
[psi,theta,phi] = EulerAngles(X_local,Y_local);
psi(isnan(psi)) = 0;
phi(isnan(phi)) = pi/2;
psi = real(psi);
% ========================= Matrizes de rotação ========================= %
M_rot = eul2rotm( [psi,theta,phi] , 'ZYX' );
M_rot = permute(M_rot,[2 1 3]);

%% Outputs %%
elem.conect = conectividade;

end

%% Função EulerAngles() %%
function [psi,theta,phi] = EulerAngles(X,Y)
% alpha = acos(-Z(:,2)./sqrt(1-Z(:,3).^2));
% beta = acos(Z(:,3));
% gama = acos(Y(:,3)./sqrt(1-Z(:,3).^2));
psi = asin(X(:,2)./sqrt(1-X(:,3).^2));
theta = asin(-X(:,3));
phi = asin(Y(:,3)./sqrt(1-X(:,3).^2));
end

%% Função RotationMatrix_ZYX() %%
function [M_rot] = RotationMatrix_ZYX(alpha,beta,gama)
s1 = sin(alpha);
s2 = sin(beta);
s3 = sin(gama);
c1 = cos(alpha);
c2 = cos(beta);
c3 = cos(gama);

s1 = reshape(s1,1,1,[]);
s2 = reshape(s2,1,1,[]);
s3 = reshape(s3,1,1,[]);
c1 = reshape(c1,1,1,[]);
c2 = reshape(c2,1,1,[]);
c3 = reshape(c3,1,1,[]);

M_rot = [ c1.*c3-c2.*s1.*s3   -c1.*s3-c2.*c3.*s1   s1.*s2;
          c3.*s1+c1.*c2.*s3   c1.*c2.*c3-s1.*s3    -c1.*s2;
          s2.*s3              c3.*s2               c2     ];


end