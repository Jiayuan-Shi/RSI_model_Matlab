function [rec_p] = self(a, b, c)
    % 矩阵维度
    n = a*(2^b) + b*2^(b-1);
    % 生成 T0 种群单倍型频率
    A0 = [ones(a, 1)/a; zeros(n-a, 1)];
    % 迭代次数
    times = c;
    % 生成 T0 个体基因型矩阵
    ind0 = (A0*A0').*~eye(n, n);
    ind0 = ind0/sum(sum(ind0));
    %% 构造总群体亲和性矩阵
    % 构造旧群体亲和矩阵
    CI_M = repmat(~eye(a, a), 2^b, 2^b);
    % 构造花粉授粉新 S-RNase 矩阵
    O1_M = zeros(2^b, b);
    for i = 1:2^b
        O1_M(i, :) = bitget(i-1, 1:b);
    end
    SI_M = kron(O1_M, ones(a, 2^(b-1)));
    % 构造新群体内部的杂交矩阵
    if b == 1
        inn_M = [0];
        O2_M = [0];
    else
        test_M = zeros(2^(b-1), (b-1));
        O2_M = zeros(2^(b-1), (b-1));
        for i = 1:2^(b-1)
            test_M(i, :) = bitget(i-1, 1:(b-1));
        end
        for i = 1:b
            O2_M((i-1)*2^(b-1)+1:i*2^(b-1), find(1:b ~= i)) = test_M;
        end
        O2_M = kron(~eye(b, b), ones(2^(b-1), 1)).*O2_M;
        inn_M = kron(O2_M, ones(1, 2^(b-1)));
    end
    % 构造新花粉授粉新 S-RNase 矩阵
    new_M = ones(b*2^(b-1), a*2^b);
    % 合并
    CI_SI_M = [CI_M, SI_M; new_M, inn_M];
    %% 构造 SLF 矩阵
    O1_exp_M = kron(O1_M, ones(a, 1));
    SLF_M = [O1_exp_M; O2_M];
    % 初始化单倍型频率输出
    rec_p = zeros(1, n);
    p = zeros(n, 1);
    % 初始化递归调用 F1
    F0 = ind0;
    for t = 1:times
        F1_i = F0;
        [F1_o, p] = F0_2_F1_ind(a, b, F1_i, CI_SI_M, SLF_M, O1_M, O2_M);
        rec_p = [rec_p; p];
        F0 = F1_o;
    end
    rec_p = rec_p(2:end, :)';
end
