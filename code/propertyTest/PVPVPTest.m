clear;
clc;

n=5;
A=rand(n);
A=abs(A);
B=diag(rand([n,1])-0.5*ones(n,1));

m=5;

sum = zeros(n,n);
for q=1:m
    for r=q:m
        sum=sum+(A^q)*B*(A^(r-q)*B*A^(m-r));
    end
end

count=0;
for i =1:n
    for j=1:n
        if sum(i,j)<(-1*10^(-5))
            disp("bad");
            count=count+1;
        end
    end
end
disp("-------")
disp(count)
disp("--end--")
        