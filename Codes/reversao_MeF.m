%% rotinas de aeroelasticidade
% A ideia é gerar uma rotina que reitere as cargas na longarina(ou
% estrutura) continuamente, para verificar a deformação e a eficiência de
% uma superfície de comando
% 
% a ideia é usar um conjunto de funções pra fazer o seguinte:
% 1- gerar uma longarina
% 2- calcular as cargas iniciasi nessa longarina (simplificado strip
% theory, talvez algo mais robusto mais pra frente)
% 3- aplicar as cargas e verificar as deformações
% 4- usar as deformações para gerar novas cargas
% 5- aplicar as novas cargas na estrutura(já deformada, alguns erros podem
% existir aqui, não sei como funciona o MEF Fontas) 
% obs: colocar um epsilon não muito pequeno e limitar o numero de iterações

%  Lista de Funções
%  ** forcas : calcula as forças aplicadas em diversas seções e suas
%  posições (entradas: coeficientes de sustentação) (saída: vetor de forças
%  e pontos)
%
%  ** ETT_FEM : calcula os deslocamentos da estrutura (entradas: longarina,
%  vetor de forças e pontos) (saída: vetor de deslocamentos(entre outros)
%
%  ** main' : realiza a reitaração das duas funções anterirores (entradas:
%  coeficientes de sutentação, longarina)  (saída: deslocamentos finais,
%  eficiência)
%
%  coef = struct com parâmetros: .cl_zero  .cl_alfa
%  aileron = struct com parâmetros: .cl_delta(cl gerado pelo aileron defletido)  
%  .pct(porcentagem da corda que corresponde ao aileron  
%   .begin e .end (ponto onde começa e onde termina o aileron_)
%   .area = area do aileron
%
%
%
%
function [eficiencia] = reversao_MeF(coef, aileron, n_pontos, b, c, v, alfa, S, L_w)
pontos = {};
forcas = {};
cl_zero = coef.cl_zero;
cl_alfa = coef.cl_alfa;
cl_delta = aileron.cl_delta;
rho = 1 ;%densidade do ar 
cargas = struct;
for i = 2: n_pontos
    pontos{i, 2} = b*(i+(i-1))/(n_pontos*2);
    pontos{i, 3} = 0.2;
    pontos{i, 1} = 0.25*c;
    forcas{i, 1} = 0;
    forcas{i, 2} = 0;
    forcas{i, 3} = (cl_zero + cl_alfa*alfa)*rho*S*v*v/2;
end

pontos{n_pontos+1, 1} = (1-aileron.pct)*c;
pontos{n_pontos+1, 2} = aileron.begin;
forcas{n_pontos+1, 3} = cl_delta*rho*aileron.area*v*v/4;
pontos{n_pontos+2, 1} = (1-aileron.pct)*c;
pontos{n_pontos+2, 2} = aileron.end;
forcas{n_pontos+2, 3} = cl_delta*rho*aileron.area*v*v/4;


cargas.F = forcas;
cargas.P = pontos;
ang_atk = zeros(1, n_pontos);
L_w = ETT_Main_FEM( L_w , cargas );

for i = 1: 50 %numero aleatorio de iterações
    for j = 1: n_pontos
        ang_atk(j) = ang_atk(j)+L_w.u(2,j, :);
        forcas{j,3} = (cl_zero + cl_alfa*(alfa+ang_atk(j))*rho*S*v*v/2);
        L_w = ETT_Main_FEM( L_w , cargas );
    end
end

eficiencia = (cl_zero + cl_alfa*(alfa+ang_atk(i))*rho*S*v*v/2)/((cl_zero + cl_alfa*alfa)*rho*S*v*v/2);

end