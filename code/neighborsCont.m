function nblist = neighborsCont(list, size, cur_nd)
    count_nb = 0;
    nblist = zeros(size,1);
    for i = 1:size
        if list(i,1)==cur_nd% || list(i,2)==cur_n
            count_nb=count_nb+1;
            nblist(count_nb) =list(i,2);
        elseif list(i,2)==cur_nd
            count_nb=count_nb+1;
            nblist(count_nb) =list(i,1);
        end
    end
    nblist = unique(nblist(1:count_nb,:));
end