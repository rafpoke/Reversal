coef = struct;
coef.cl_zero = 1;
coef.cl_alfa = 1;

aileron = struct;
aileron.cl_delta = 1;
aileron.pct = 1;
aileron.area = 1;
aileron.begin = 1;
aileron.end = 1;

fileName = 'longarina_full.mat'; %longarina que vai ser usada
load(fileName);

Longarina = L_w;

b = 1;
c = 1;
v = 1;
alfa = 1;
S = 1;
n_pontos = 200;

eficiencia = reversao_MeF(coef, aileron, n_pontos, b, c, v, alfa, S, Longarina);