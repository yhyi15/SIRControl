function inCont = findedgeCont(list, size_l, snd, tnd)
    inCont = false;
    for it = 1:size_l
        if ((snd==list(it,1) && tnd == list(it,2)) || (snd==list(it,2) && tnd == list(it,1)))
            inCont = true;
            break;
        end
    end
end