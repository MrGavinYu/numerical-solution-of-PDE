function [heat, x, t] = CN(rho, c, k, height, tao, time_span, h, Tin, Tout, Tinit)
% CN ����ԳƸ�ʽ���
% Gavin <www.bigbugs.cn>

% ���ֿռ�����
[x, area] = DevideX(height, h);
t = 0:tao:time_span;
a = k./(rho.*c);
% �ҳ��ֽ���
L = cumsum(height);
bound = zeros(1, length(height)-1);
for i = 1:length(bound)
    bound(i) = find(L(i) <= x+eps & L(i) >= x-eps);
end
% ��ÿ���������֮��Ĳ���
h = diff(x);
h = [h(1) h];
heat = zeros(length(t), length(x));
% ��ֵ����
heat(1, :) = Tinit;
for m = 2:length(t)
    K = sparse(length(x), length(x));
    b = zeros(length(x), 1);
    for j = 2:length(x)-1 % �Ե�j��x�������
        M = 2/(h(j+1) + h(j))/h(j+1);
        N = 2/(h(j+1) + h(j))/h(j);
        pos = area(j);
        K(j, j-1) = -a(pos)*tao*N/2;
        K(j, j) = 1 + a(pos)*tao*M/2 + a(pos)*tao*N/2;
        K(j, j+1) = -a(pos)*tao*M/2;
        b(j) = a(pos)*tao*N/2 * heat(m-1, j-1) +  (a(pos)*tao*(-M-N)/2 + 1) * heat(m-1, j)  + a(pos)*tao*M/2 * heat(m-1, j+1); 
    end
    % ���ӱ�ֵ
    K(1, 1) = 1;K(end, end) = 1;
    b(1) = Tin;b(end) = Tout;
    % ���ʽ����������
    for i = 1:length(bound)
        j = bound(i);
        pos = area(j);
        K(j, j-1) = -k(pos)/h(j);
        K(j, j) = k(pos)/h(j) + k(pos+1)/h(j+1);
        K(j, j+1) = -k(pos+1)/h(j+1);
        b(j) = 0;
    end
    heat(m, :) = K\b;
end
end

