function [ Longarina ] = ETT_PontosTirante( geom , Longarina )
Malha = {geom.Malha};
is_long = [Longarina.is_tirante];
% fcC = {Longarina.fcC};
fc = {Longarina.fc};
fb = {Longarina.fb};
% fhC = {Longarina.fhC};
num_pontos = [Longarina.num_pontos];
idx_malha = {Longarina.idx_malha};

num_long = length(Longarina);

%% CORPO DA FUNÇÃO %%
aux = 1:num_long;
aux = aux(is_long);

% Para simplificar os programas, cada longarina é tratada separadamente
for i = aux
    [Longarina(i).pontos , Longarina(i).posLong] = pontosTiranteSimples( {Malha{idx_malha{i}}} , fb{i} , fc{i} , num_pontos(i) );
end
end


%% Função pontosTiranteSimples() %%
function [ pontos_longarina , posLong ] = pontosTiranteSimples( Malha , fb , fc , num_pontos )
%% Carregando Valores %%
Malha1 = Malha{1};
Malha2 = Malha{2};

X1 = Malha1.X;
Y1 = Malha1.Y;
Z_med1 = Malha1.Z_med;

X2 = Malha2.X;
Y2 = Malha2.Y;
Z_med2 = Malha2.Z_med;

b = Y1(1,end);                                                               % semi-envergadura da malha

%% Corpo da Função %%
% ============================ INICIALIZAÇÃO ============================ %
pontos_longarina = zeros(3,num_pontos);
y_long_disc = zeros(1,2);

% ======================== Posição em Y discreta ======================== %
y_long_disc = fb * b;

% ====================== Posição em Y da longarina ====================== %
y_long = linspace(y_long_disc(1),y_long_disc(end),num_pontos);
pontos_longarina(2,:) = y_long;

% =============================== posLong =============================== %
posLong = [1 num_pontos];

% =================== CORDA NAS POSIÇÕES DA LONGARINA =================== %
cordas = [ interp1(Y1(1,:),X1(end,:)-X1(1,:),y_long(1)) , interp1(Y2(1,:),X2(end,:)-X2(1,:),y_long(2)) ];

% ====================== Posição em X da longarina ====================== %
dx = fc .* cordas;
x0 = [ interp1(Y1(1,:),X1(1,:),y_long(1)) , interp1(Y2(1,:),X2(1,:),y_long(2)) ];
x_long = linspace(x0(1)+dx(1),x0(2)+dx(2),num_pontos);
pontos_longarina(1,:) = x_long;

% ====================== Posição em Z da longarina ====================== %
z = [ interp1(Y1(1,:),Z_med1(1,:),y_long(1)) , interp1(Y2(1,:),Z_med2(1,:),y_long(2)) ];
z_long = linspace(z(1),z(2),num_pontos);
pontos_longarina(3,:) = z_long;

pontos_longarina = pontos_longarina';

end
