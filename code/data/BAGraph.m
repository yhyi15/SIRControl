%% generates a power law graph using the
% Barabasi-Albert (BA) model
% http://en.wikipedia.org/wiki/BA_Model

function A = BAGraph(N,m,m0)


% each node has weighted probability of getting new edge
% deg(i) / sum_j (deg(j))

% start with empty adjaceny matrix
A = zeros(N,N);



% start with first m0 nodes connected in a line
for ndx=2:m0
    A(ndx-1,ndx) = 1;
    A(ndx,ndx-1) = 1;
end

% keep track of node degrees
deg = sum(A);

% add the remaining nodes one at a time.
for i=m0+1:N
    disp(i);
    % add m edges per new node
    j=1;
    while j<=m
        r = rand;
        denom = sum(deg);
        prob = deg ./ denom;
        p = 0;
        for k=1:i-1
            p = p + prob(k);
            if (r <= p && A(i,k)==0)
                % connect new node(i) to node k
                A(i,k) = 1;
                A(k,i) = 1;
                deg = sum(A);
                j=j+1;
                break;
            end
            if (r<=p && A(i,k)==1)
                break;
            end
        end
    end
end
end
