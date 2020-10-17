%takes rows and columns of the current inverse, returns marginal gain
function [val]=numUpdate(nn, Ri,Rj,Ci,Cj, x0, r0, lft, rt, ei, ej, Aij, Aji)
    cij = (1-x0(ei)-r0(ei))*Aij;
    LUpdateij = Ci;%curInv(:, ei);
    RUpdateij = Rj;%curInv(ej, :);
    numeratorij = 1 + cij * Rj(ei);%curInv(ej,ei);
    val1 = lft*rt- cij * 1/numeratorij * (lft*LUpdateij) * (RUpdateij*x0);
    cji = (1-x0(ej)-r0(ej))*Aji;
    LUpdateji = Cj -  cij * 1/numeratorij * LUpdateij * RUpdateij(ej);
    RUpdateji = Ri -  cij * 1/numeratorij * RUpdateij * LUpdateij(ei);
    updatedji = Ri(ej) - cij * 1/numeratorij * LUpdateij(ei) * RUpdateij(ej);
    numeratorji = 1 + cji * updatedji;%curInv(ei,ej);
    val = val1 - cji * 1/numeratorji * (lft*LUpdateji) * (RUpdateji*x0);
    
end