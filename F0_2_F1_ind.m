function [F1, p] = F0_2_F1_ind(a, b, F0, CI_SI_M, SLF_M, O1_M, O2_M)
    % 群体大小
    popsize = 5000;
    % 变异概率
    pm = 0.05;
    pe = 0.05;
    % 矩阵维度
    n = a*(2^b) + b*2^(b-1);
    % 初始化单倍型频率 A
    A = sum(F0)';
    % 初始化 index
    index = 1:n;
    index(a*2^b:2^(b-1):n) = [];
    % 突变及重组影响到生殖事件的花粉单倍型频率
    for m = index
        for o = index
            % 累积新的单倍型频率作为授粉时的输入
            acc_m = find_slf(a, b, O1_M, O2_M, m);
            acc_o = find_slf(a, b, O1_M, O2_M, o);
            accum = [acc_m, acc_o];
            % 累积
            if length(accum) ~= 0
                A(accum) = A(accum) + (F0(m, o))*binornd(1, pe/(length(accum)), length(accum), 1)/(2*popsize);
            end
        end
    end
    A = A/sum(A);
    % 突变 S-RNase
    % 获得 SLF_M 和 O2_M 行状态
    SLF_ind = cellfun(@(row) num2str(row, '%d'), num2cell(SLF_M, 2), 'UniformOutput', false);
    SLF_ind = cellfun(@bin2dec, SLF_ind);
    O2_ind = cellfun(@(row) num2str(row, '%d'), num2cell(O2_M, 2), 'UniformOutput', false);
    O2_ind = cellfun(@bin2dec, O2_ind);
    % 初始化
    A0 = A;
    for d = 0:2^b-1
        co_SLF = sum(A(find(d == SLF_ind)));
        subjct = find(d == O2_ind);
        if length(subjct) ~= 0
            A0(a*2^b+subjct) = A0(a*2^b+subjct) + co_SLF*(binornd(1, pm/length(subjct), length(subjct), 1))/(2*popsize);
        end
    end
    A = A0/sum(A0);
    % 初始化 Fn+1 子代的基因型矩阵
    fi = zeros(size(F0));
    % 遍历母本基因型为 Si/Sj
    for i = 1:n
        for j = 1:n
            % 模仿授粉事件进行随机扰动
            tmp_A = A.*(1+0.05*randn(n,1)');
            % 等比例扩增 tmp_A 使和为 1, 作为输入
            tmp_A = tmp_A/sum(tmp_A);
            % 单次授粉事件的花粉单倍型频率计算
            acc = CI_SI_M(:,i) & CI_SI_M(:,j);
            x = acc.*tmp_A;
            if sum(x) ~= 0
                x = x/sum(x);
            end
            % 获得 Fn+1 个体基因型矩阵,i,j 基因型分别和花粉单倍型组合
            fi(i, :) = fi(i, :) + F0(i, j).*x';
            fi(j, :) = fi(j, :) + F0(i, j).*x';
        end
    end
    % 基因型矩阵均一化处理
    fi = fi + fi';
    fi = fi/sum(sum(fi));
    % 返回子代基因型频率 F1
    F1 = fi;
    % 单倍型频率均一化处理
    A = sum(fi);
    A = A/sum(A);
    % 返回子代单倍型频率 A
    p = A;
end
