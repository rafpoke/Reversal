%% rotinas de aeroelasticidade
% A ideia � gerar uma rotina que reitere as cargas na longarina(ou
% estrutura) continuamente, para verificar a deforma��o e a efici�ncia de
% uma superf�cie de comando
% 
% a ideia � usar um conjunto de fun��es pra fazer o seguinte:
% 1- gerar uma longarina
% 2- calcular as cargas iniciasi nessa longarina (simplificado strip
% theory, talvez algo mais robusto mais pra frente)
% 3- aplicar as cargas e verificar as deforma��es
% 4- usar as deforma��es para gerar novas cargas
% 5- aplicar as novas cargas na estrutura(j� deformada, alguns erros podem
% existir aqui, n�o sei como funciona o MEF Fontas) 
% obs: colocar um epsilon n�o muito pequeno e limitar o numero de itera��es

%  Lista de Fun��es
%  ** forcas : calcula as for�as aplicadas em diversas se��es e suas
%  posi��es (entradas: coeficientes de sustenta��o) (sa�da: vetor de for�as
%  e pontos)
%
%  ** ETT_FEM : calcula os deslocamentos da estrutura (entradas: longarina,
%  vetor de for�as e pontos) (sa�da: vetor de deslocamentos(entre outros)
%
%  ** main' : realiza a reitara��o das duas fun��es anterirores (entradas:
%  coeficientes de sutenta��o, longarina)  (sa�da: deslocamentos finais,
%  efici�ncia)
%
%  coef = struct com par�metros: .cl_zero  .cl_alfa
%  aileron = struct com par�metros: .cl_delta(cl gerado pelo aileron defletido)  
%  .pct(porcentagem da corda que corresponde ao aileron  
%   .begin e .end (ponto onde come�a e onde termina o aileron_)
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

for i = 1: 50 %numero aleatorio de itera��es
    for j = 1: n_pontos
        ang_atk(j) = ang_atk(j)+L_w.u(2,j, :);
        forcas{j,3} = (cl_zero + cl_alfa*(alfa+ang_atk(j))*rho*S*v*v/2);
        L_w = ETT_Main_FEM( L_w , cargas );
    end
end

eficiencia = (cl_zero + cl_alfa*(alfa+ang_atk(i))*rho*S*v*v/2)/((cl_zero + cl_alfa*alfa)*rho*S*v*v/2);

end