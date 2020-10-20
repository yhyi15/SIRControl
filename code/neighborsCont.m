function nblist = neighborsCont(list, size, cur_nd)
    count_nb = 0;
    nblist = zeros(size, 2);
    for i = 1:size
        if list(i,1)==cur_nd || list(i,2)==cur_nd
            count_nb=count_nb+1;
            nblist(count_nb,1) =list(i,1);
            nblist(count_nb,2) =list(i,2);
        end
    end
    nblist = nblist(1:count_nb,:);
end