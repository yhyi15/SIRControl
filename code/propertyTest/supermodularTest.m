clear;
clc;

n = 30;
delta = 0.25;

for round=1:100
    A=zeros(n,n);
    %A1=A; A2=A; Aint=A; Auni=A;
    G1 = graph(A);
    G2 = graph(A);
    Gint = graph(A);
    Guni = graph(A);
    for i = 1:n-1
        for j = i+1:n
            edgeCoin1 = binornd(1,0.3);
            edgeCoin2 = binornd(1,0.2);
            if edgeCoin1 ==1
                G1 = addedge(G1, i,j, 1);
            end
            if edgeCoin2 ==1
                G2 = addedge(G2, i,j, 1);
            end
            if edgeCoin1 ==1 && edgeCoin2==1
                Gint = addedge(Gint, i,j, 1);
            end
            if edgeCoin1 ==1 || edgeCoin2 ==1
                Guni = addedge(Guni, i, j,1);
            end
        end
    end
    A1 = adjacency(G1,'weighted');
    A2 = adjacency(G2,'weighted');
    Aint = adjacency(Gint,'weighted');
    Auni = adjacency(Guni,'weighted');
    
    
    dmax1 = max(sum(A1));
    dmax2 = max(sum(A2));
    
    
    dmax = max(dmax1+dmax2);
    
    beta = delta/(dmax+1);
    
    betaList = beta* ones(n,1);
    deltaList = delta * ones(n,1);
    B = diag(betaList);
    D = diag(deltaList);
    
    x0 = zeros(n, 1);
    r0 = zeros(n, 1);
    
    %set initial conditions
    s = 5;
    %choose s seeds
    S = randsample(n,s);
    % initiate x and r
    for i = 1: s
        x0(i) = rand()/3;
    end
    for i = 1:n
        r0(i) = rand()/100;
    end
    
    X0 = diag(x0);
    R0 = diag(r0);
    I = eye(n);
    M1 = I - D + (I-X0-R0)*B*A1;
    M2 = I - D + (I-X0-R0)*B*A2;
    Mint = I - D + (I-X0-R0)*B*Aint;
    Muni = I - D + (I-X0-R0)*B*Auni;
    
    
    sigma1 = ones(1,n)* (M1+D-I) * ((I-M1)\x0);
    sigma2 = ones(1,n)* (M2+D-I) * ((I-M2)\x0);
    sigmauni = ones(1,n)* (Muni/2+D-I) * ((I-Muni/2)\x0);
    sigmaint = ones(1,n)* (Mint/2+D-I) * ((I-Mint/2)\x0);
    
    if (sigma1+sigma2)<(sigmaint+sigmauni-10^-6)
        disp("nonconvex");
    end
    disp(round);
    disp("end");
end