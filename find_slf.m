function [f2] = find_slf(a, b, O1_M, O2_M, ij)
    % Determine operation based on the range of ij
    if ij < a*2^b
        sel = fix((ij-1)/a) + 1;
        sum_rows = sum(O1_M');
        sum_sel = sum(O1_M(sel,:));
        test1 = [];
        rows = 2^b;
        for m = 1:rows
            % Check for same existing genotype and new SLF count is one more than before
            test1 = [test1, (isequal((O1_M(sel,:) & O1_M(m,:)), O1_M(sel,:))) & (sum_rows(m) == sum_sel + 1)];
        end
        ref = find(test1 == 1) - 1;
        re = mod(ij, a);
        if re == 0
            re = a;
        end
        f2 = a*ref + re;
    elseif ij >= a*2^b + 1 & ij < a*2^b + b*2^(b-1)
        sel = ij - a*2^b;
        blo = fix((sel-1)/2^(b-1));
        re = mod(sel, 2^(b-1));
        if re == 0
            re = 2^(b-1);
        end
        sel_M = O2_M(blo*2^(b-1)+1 : blo*2^(b-1)+2^(b-1), :);
        sum_rows = sum(sel_M');
        sum_sel = sum(sel_M(re, :));
        test1 = [];
        rows = 2^(b-1);
        for m = 1:rows
            test1 = [test1, (isequal((sel_M(re,:) & sel_M(m,:)), sel_M(re,:))) & (sum_rows(m) == sum_sel + 1)];
        end
        ref = find(test1 == 1);
        f2 = ref + a*2^b + blo*2^(b-1);
    else
        f2 = 0;
    end
end
