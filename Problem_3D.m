% Основные параметры задачи
N = 10;
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

tic;

% задание точного xi - exactXi
exactXi = zeros(N, N, N);
xi = zeros(N^3, 1);

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
  
% задаём матрицы a
% счет индексов метода квадратур
ind = zeros(1, N);
for k = 1:N
   ind(k) = 1.0;
end
ind(1) = 0.5;
ind(N) = 0.5;
    
a_10 = zeros(N^3); 
a_20 = zeros(N^3); 
a_30 = zeros(N^3); 
for k = 1:N
    for l = 1:N
        for m = 1:N
            coord_1 = N^2 * (k - 1) + N * (l - 1) + m;
            for p = 1:N
                for q = 1:N
                    for r = 1:N
                        coord_2 = N^2 * (p - 1) + N * (q - 1) + r;
                        a_10(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([h * (k - 1) h * (l - 1) h * (m - 1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                        a_20(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([h * (k - 1) h * (l - 1) h * (m - 1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                        a_30(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([h * (k - 1) h * (l - 1) h * (m - 1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                    end
                end
            end
        end
    end
end

% вычисляем матрицу over_a для каждого источника
a_11 = zeros(N^3);
a_12 = zeros(N^3);
a_13 = zeros(N^3);
a_14 = zeros(N^3);
a_15 = zeros(N^3);
a_16 = zeros(N^3);
a_21 = zeros(N^3);
a_22 = zeros(N^3);
a_23 = zeros(N^3);
a_24 = zeros(N^3);
a_25 = zeros(N^3);
a_26 = zeros(N^3);
a_31 = zeros(N^3);
a_32 = zeros(N^3);
a_33 = zeros(N^3);
a_34 = zeros(N^3);
a_35 = zeros(N^3);
a_36 = zeros(N^3);

 for k = 1:N
     for l = 1:N
         angle_S_phi_1 = Detect_Angle_phi(1, 1) + h_ang_phi * (l - 1);
         angle_S_phi_2 = Detect_Angle_phi(2, 1) + h_ang_phi * (l - 1);
         angle_S_phi_3 = Detect_Angle_phi(3, 1) + h_ang_phi * (l - 1);
         angle_S_phi_4 = Detect_Angle_phi(4, 1) + h_ang_phi * (l - 1);
         angle_S_phi_5 = Detect_Angle_phi(5, 1) + h_ang_phi * (l - 1);
         angle_S_phi_6 = Detect_Angle_phi(6, 1) + h_ang_phi * (l - 1);
         for m = 1:N
             coord_1 = N^2 * (k - 1) + N * (l - 1) + m;
             angle_S_thetha_1 = Detect_Angle_thetha(1, 1)  + h_ang_thetha * m;
             angle_S_thetha_2 = Detect_Angle_thetha(2, 1)  + h_ang_thetha * m;
             angle_S_thetha_3 = Detect_Angle_thetha(3, 1)  + h_ang_thetha * m;
             angle_S_thetha_4 = Detect_Angle_thetha(4, 1)  + h_ang_thetha * m;
             angle_S_thetha_5 = Detect_Angle_thetha(5, 1)  + h_ang_thetha * m;
             angle_S_thetha_6 = Detect_Angle_thetha(6, 1)  + h_ang_thetha * m;
             for p = 1:N
                 for q = 1:N
                     for r = 1:N  
                         coord_2 = N^2 * (p - 1) + N * (q - 1) + r;
                         a_11(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * cos(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * sin(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_21(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * cos(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * sin(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_31(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * cos(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_1) * sin(angle_S_phi_1) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_1)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                         a_12(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * cos(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * sin(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_2)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_22(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * cos(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * sin(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_2)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_32(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * cos(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_2) * sin(angle_S_phi_2) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_2)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                         a_13(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * cos(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * sin(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_3)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_23(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * cos(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * sin(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_3)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_33(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * cos(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_3) * sin(angle_S_phi_3) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_3)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                         a_14(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * cos(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * sin(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_4)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_24(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * cos(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * sin(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_4)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_34(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * cos(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_4) * sin(angle_S_phi_4) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_4)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                         a_15(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * cos(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * sin(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_5)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_25(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * cos(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * sin(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_5)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_35(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * cos(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_5) * sin(angle_S_phi_5) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_5)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                         a_16(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * cos(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * sin(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_6)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(1), c_0);
                         a_26(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * cos(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * sin(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_6)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(2), c_0);
                         a_36(coord_1, coord_2) = h^3 * ind(p) * ind(q) * ind(r) * G([0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * cos(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * sin(angle_S_thetha_6) * sin(angle_S_phi_6) 0.5 + (R_sourse + h_D * (k - 1)) * cos(angle_S_thetha_6)], [h * (p - 1) h * (q - 1) h * (r - 1)], omega(3), c_0);
                     end
                 end
             end
         end
     end
 end

% для правой части основной системы
f_1_R_1 = zeros(N^3, 1);
f_1_R_2 = zeros(N^3, 1);
f_1_R_3 = zeros(N^3, 1);
f_1_R_4 = zeros(N^3, 1);
f_1_R_5 = zeros(N^3, 1);
f_1_R_6 = zeros(N^3, 1);

f_2_R_1 = zeros(N^3, 1);
f_2_R_2 = zeros(N^3, 1);
f_2_R_3 = zeros(N^3, 1);
f_2_R_4 = zeros(N^3, 1);
f_2_R_5 = zeros(N^3, 1);
f_2_R_6 = zeros(N^3, 1);

f_3_R_1 = zeros(N^3, 1);
f_3_R_2 = zeros(N^3, 1);
f_3_R_3 = zeros(N^3, 1);
f_3_R_4 = zeros(N^3, 1);
f_3_R_5 = zeros(N^3, 1);
f_3_R_6 = zeros(N^3, 1);

for k = 1:N
    for l = 1:N
        for m = 1:N
            coord_1 = N^2 * (k -1) + N * (l - 1) + m;
            f_1_R_1(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(1, :), omega(1), c_0);
            f_1_R_2(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(2, :), omega(1), c_0);
            f_1_R_3(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(3, :), omega(1), c_0);
            f_1_R_4(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(4, :), omega(1), c_0);
            f_1_R_5(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(5, :), omega(1), c_0);
            f_1_R_6(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(6, :), omega(1), c_0);
            
            f_2_R_1(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(1, :), omega(2), c_0);
            f_2_R_2(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(2, :), omega(2), c_0);
            f_2_R_3(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(3, :), omega(2), c_0);
            f_2_R_4(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(4, :), omega(2), c_0);
            f_2_R_5(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(5, :), omega(2), c_0);
            f_2_R_6(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(6, :), omega(2), c_0);
            
            f_3_R_1(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(1, :), omega(3), c_0);
            f_3_R_2(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(2, :), omega(3), c_0);
            f_3_R_3(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(3, :), omega(3), c_0);
            f_3_R_4(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(4, :), omega(3), c_0);
            f_3_R_5(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(5, :), omega(3), c_0);
            f_3_R_6(coord_1) = G([h * (k - 1) h * (l - 1) h * (m - 1)], Source(6, :), omega(3), c_0);
            
            
        end
    end
end

% Вспомогательная матрица для первой частоты
Matr = zeros(N^3);
for ii = 1:N^3
    for jj = 1:N^3
        Matr(ii, jj) = -a_10(ii, jj) * xi(jj);
    end
    Matr(ii, ii) = Matr(ii, ii) + 1;
end
u11 = Matr \ f_1_R_1;
f_1_S_1 = a_11 * (xi .* u11);
u12 = Matr \ f_1_R_2;
f_1_S_2 = a_12 * (xi .* u12);
u13 = Matr \ f_1_R_3;
f_1_S_3 = a_13 * (xi .* u13);
u14 = Matr \ f_1_R_4;
f_1_S_4 = a_14 * (xi .* u14);
u15 = Matr \ f_1_R_5;
f_1_S_5 = a_15 * (xi .* u15);
u16 = Matr \ f_1_R_6;
f_1_S_6 = a_16 * (xi .* u16);

% Вспомогательная матрица для второй частоты
for ii = 1:N^3
    for jj = 1:N^3
        Matr(ii, jj) = -a_20(ii, jj) * xi(jj);
    end
    Matr(ii, ii) = Matr(ii, ii) + 1;
end
u21 = Matr \ f_2_R_1;
f_2_S_1 = a_21 * (xi .* u21);
u22 = Matr \ f_2_R_2;
f_2_S_2 = a_22 * (xi .* u22);
u23 = Matr \ f_2_R_3;
f_2_S_3 = a_23 * (xi .* u23);
u24 = Matr \ f_2_R_4;
f_2_S_4 = a_24 * (xi .* u24);
u25 = Matr \ f_2_R_5;
f_2_S_5 = a_25 * (xi .* u25);
u26 = Matr \ f_2_R_6;
f_2_S_6 = a_26 * (xi .* u26);

% Вспомогательная матрица для третьей частоты
for ii = 1:N^3
    for jj = 1:N^3
        Matr(ii, jj) = -a_30(ii, jj) * xi(jj);
    end
    Matr(ii, ii) = Matr(ii, ii) + 1;
end
u31 = Matr \ f_3_R_1;
f_3_S_1 = a_31 * (xi .* u31);
u32 = Matr \ f_3_R_2;
f_3_S_2 = a_32 * (xi .* u32);
u33 = Matr \ f_3_R_3;
f_3_S_3 = a_33 * (xi .* u33);
u34 = Matr \ f_3_R_4;
f_3_S_4 = a_34 * (xi .* u34);
u35 = Matr \ f_3_R_5;
f_3_S_5 = a_35 * (xi .* u35);
u36 = Matr \ f_3_R_6;
f_3_S_6 = a_36 * (xi .* u36);

%---------------------------------------
toc;
%---------------------------------------