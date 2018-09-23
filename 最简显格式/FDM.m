function [heat, x, t] = FDM(rho, c, k, height, tao, time_span, h, Tin, Tout, Tinit)
% FDM �ȴ���������ǰ���
% Gavin <www.bigbugs.cn>

L = cumsum(height);
a = k./(rho.*c);
% ����ȶ�������
r = a.*tao./(h.^2);
if(sum(r > 0.5))
    error('��%d�㲽�����ò������ȶ���Ҫ��:�����=%f\n', find(r>0.5), r(r>0.5))
    % ������ִ˴���Ҫ����Сʱ�䲽����Ҫ������ռ䲽����
    % ʹ�������С��0.5�����ȶ�������
end
% ���ֿռ�����
[x, area] = DevideX(height, h);
t = 0:tao:time_span;
% ��ÿ���������֮��Ĳ���
h = diff(x);
h = [h(1) h];
% �ҳ����Ϸֽ���
bound = zeros(1, length(height)-1);
for i = 1:length(bound)
    bound(i) = find(L(i) <= x+eps & L(i) >= x-eps);
end
% ��ֵ����
heat = zeros(length(t), length(x));
heat(1, :) = Tinit;
for i = 2:length(t)
    heat(i, 1) = Tin;% ��һ���һ���ֵ����
    heat(i, end) = Tout;% ���һ���һ���ֵ����
    for j = 2:length(x)-1
        heat(i, j) = (heat(i-1, j+1)-heat(i-1, j))/h(j+1)/(h(j+1) + h(j)) - (heat(i-1, j)-heat(i-1, j-1))/h(j)/(h(j+1)+h(j));
        heat(i, j) = heat(i, j)*k(area(j))/rho(area(j))/c(area(j)) * tao + heat(i-1, j);
    end
    % �ֽ�������
    for j = 1:length(bound)
        b = bound(j);
        k_l = k(area(b-1));
        k_n = k(area(b+1));
        a = k_l/k_n*h(b+1)/h(b);
        heat(i, b) = (a*heat(i, b-1) + heat(i, b+1))/(a+1);
    end
end
end