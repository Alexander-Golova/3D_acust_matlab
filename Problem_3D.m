% Основные параметры задачи
N = 20;
h = 1 / N;

% параметры регуляризованного метода
numberOfIterations = 20;
alpha = 0.75;
multiplier = 0.25;

% задание характеристик поля
c_0 = 1.0;
omega = [2.0 * pi 3.0 * pi 4.0 * pi];

% Характеристики эксперимента
% удаление детектора и источников от неоднородности
rho_S = 1.0;
% ширина зоны детекторов
width_D = 1.0;
% шаг по радиусу
h_D = width_D / N;
% шаг по углу
h_ang = pi / N;

N = N + 1;

% количество источников
numberSource = 6;
% координаты источников в трехмерном пространстве
% источники расположены на сфере радиуса R_sourse
R_sourse = 0.5 * sqrt(3) + rho_S;

angle_phi = 0. * pi;
angle_thetha = 0.5 * pi;
S1 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

angle_phi = 1.0 * pi;
angle_thetha = 0.5 * pi;
S2 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

angle_phi = 0.5 * pi;
angle_thetha = 0.5 * pi;
S3 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

angle_phi = 1.5 * pi;
angle_thetha = 0.5 * pi;
S4 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

angle_phi = 0.0 * pi;
angle_thetha = 0.0 * pi;
S5 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

angle_phi = 0.0 * pi;
angle_thetha = 1.0 * pi;
S6 = [0.5 + R_sourse * sin(angle_thetha) * cos(angle_phi) 0.5 + R_sourse * sin(angle_thetha) * sin(angle_phi) 0.5 + R_sourse * cos(angle_thetha)];

Source = [S1; S2; S3; S4; S5; S6];
Source_Angle = [0.0 * pi 0.5 * pi; 1.0 * pi 0.5 * pi; 0.5 * pi 0.5 * pi; 1.5 * pi 0.5 * pi; 0.0 * pi 0.0 * pi; 0.0 * pi 1.0 * pi];
Detect_Angle_phi = [0.5 * pi 1.5 * pi; -0.5 * pi 0.5 * pi; 1.0 * pi 2.0 * pi; 0.0 * pi 1.0 * pi; 0.0 * pi 2.0 * pi; 0.0 * pi 2.0 * pi];
Detect_Angle_thetha = [0.0 * pi 1.0 * pi; 0.0 * pi 1.0 * pi; 0.0 * pi 1.0 * pi; 0.0 * pi 1.0 * pi; 0.5 * pi 1.0 * pi; 0.0 * pi 0.5 * pi];

% координаты приемников зависят от источника
% номер источника
number_Sourse = 1;
X_Res = zeros(N^3, 1);
Y_Res = zeros(N^3, 1);
Z_Res = zeros(N^3, 1);
X_pl = zeros(N^3, 1);
Y_pl = zeros(N^3, 1);
Z_pl = zeros(N^3, 1);
h_ang_phi = (Detect_Angle_phi(number_Sourse, 2) - Detect_Angle_phi(number_Sourse, 1)) / (N - 1);
h_ang_thetha = (Detect_Angle_thetha(number_Sourse, 2) - Detect_Angle_thetha(number_Sourse, 1)) / (N - 1);

u_pl = 0 : h : 1.0;
v_pl = 0 : h : 1.0;
[X_pl_2D, Y_pl_2D] = meshgrid(u_pl, v_pl);

for k = 1:N
    for l = 1:N
        for m = 1:N
            angle_S_phi = Detect_Angle_phi(number_Sourse, 1) + h_ang_phi * (l - 1);
            angle_S_thetha = Detect_Angle_thetha(number_Sourse, 1)  + h_ang_thetha * (m - 1);
            X_Res(N^2 * (k - 1) + N * (l - 1) + m) = 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha) * cos(angle_S_phi);
            Y_Res(N^2 * (k - 1) + N * (l - 1) + m) = 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha) * sin(angle_S_phi);
            Z_Res(N^2 * (k - 1) + N * (l - 1) + m) = 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha);
            X_pl(N^2 * (k - 1) + N * (l - 1) + m) = h * (k - 1);
            Y_pl(N^2 * (k - 1) + N * (l - 1) + m) = h * (l - 1);
            Z_pl(N^2 * (k - 1) + N * (l - 1) + m) = h * (m - 1);
        end
    end
end

% визуализация источников и приёмников
% создаем графическое окно
hF = figure;
% делаем нулевое графическое окно текущим
figure(hF);
plot3(X_pl, Y_pl, Z_pl, '.k', Source(number_Sourse, 1), Source(number_Sourse, 2), Source(number_Sourse, 3), 'or', X_Res, Y_Res, Z_Res, '.b')

%---------------------------------------------------------------
% 
% Решаем прямую задачу
%
%---------------------------------------------------------------

% задание точного xi - exactXi
exactXi = zeros(N, N, N);
xi = zeros(N^3, 1);
xi_pribl = zeros(N^3, 1);

for k = 1:N
    for l = 1:N
        for m = 1:N
            temp = 0.3 * exp(-((k * h - 0.7)^2 + (l * h - 0.7)^2 + (m * h - 0.5)^2) * 64);
            exactXi(k, l, m) = temp;
            xi(N^2 * (k - 1) + N * (l - 1) + m) = temp;
        end
    end
end

% рисуем точное решение
hF = figure;
figure(hF);
mesh(X_pl_2D, Y_pl_2D, exactXi(:, :, 3));
  

