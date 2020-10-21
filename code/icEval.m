function numInfect = icEval(Gr, nn, rounds,S ,s, A)
cur_count_f = 0;
for round_f = 1:10*rounds
    terminal_f = randsample(nn,1);
    
    stack_f = zeros(nn, 1);
    stack_f(1) = terminal_f;
    top_f = 1;
    activated_f = zeros(nn,1);
    
    while top_f>0
        %pop
        cur_f = stack_f(top_f);
        top_f = top_f -1;
        if activated_f(cur_f)==0
            nb_f = neighbors(Gr, cur_f);
            for iter_f = 1:size(nb_f,1)
                %push
                t_f = nb_f(iter_f);
                contCoin_f = rand();
                if (contCoin_f<=A(cur_f,t_f))
                    %push
                    top_f = top_f+1;
                    stack_f(top_f) = t_f;
                end
            end
            activated_f(cur_f)=1;
        end
    end
    for seed_iter_f= 1:s
        if activated_f(S(seed_iter_f))==1
            cur_count_f=cur_count_f+1;
            break;
        end
    end
end
numInfect = cur_count_f*nn/10/rounds;
end