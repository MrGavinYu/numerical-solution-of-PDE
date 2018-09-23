function [heat, x, t] = FDM(rho, c, k, height, tao, time_span, h, Tin, Tout, Tinit)
% FDM 热传导方程向前差分
% Gavin <www.bigbugs.cn>

L = cumsum(height);
a = k./(rho.*c);
% 检查稳定性条件
r = a.*tao./(h.^2);
if(sum(r > 0.5))
    error('第%d层步长设置不满足稳定性要求:网格比=%f\n', find(r>0.5), r(r>0.5))
    % 如果出现此错误，要不调小时间步长，要不调大空间步长，
    % 使其网格比小于0.5，即稳定性条件
end
% 划分空间网格
[x, area] = DevideX(height, h);
t = 0:tao:time_span;
% 求每两个网格点之间的步长
h = diff(x);
h = [h(1) h];
% 找出材料分界面
bound = zeros(1, length(height)-1);
for i = 1:length(bound)
    bound(i) = find(L(i) <= x+eps & L(i) >= x-eps);
end
% 初值条件
heat = zeros(length(t), length(x));
heat(1, :) = Tinit;
for i = 2:length(t)
    heat(i, 1) = Tin;% 第一层第一类边值条件
    heat(i, end) = Tout;% 最后一层第一类边值条件
    for j = 2:length(x)-1
        heat(i, j) = (heat(i-1, j+1)-heat(i-1, j))/h(j+1)/(h(j+1) + h(j)) - (heat(i-1, j)-heat(i-1, j-1))/h(j)/(h(j+1)+h(j));
        heat(i, j) = heat(i, j)*k(area(j))/rho(area(j))/c(area(j)) * tao + heat(i-1, j);
    end
    % 分界面条件
    for j = 1:length(bound)
        b = bound(j);
        k_l = k(area(b-1));
        k_n = k(area(b+1));
        a = k_l/k_n*h(b+1)/h(b);
        heat(i, b) = (a*heat(i, b-1) + heat(i, b+1))/(a+1);
    end
end
end