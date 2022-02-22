% функция Грина для задачи акустики
function y = G(P1, P2, omega, c_0)
dist = sqrt((P1(1) - P2(1))^2 + (P1(2) - P2(2))^2 + (P1(3) - P2(3))^2);
k = omega / c_0;

if (dist < 0.001)
    y = 0.0;
else
    y = exp(-1i * dist * k) / (4.0 * dist * pi);
end