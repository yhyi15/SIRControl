% update inverse
function [udtInv]=invUpdate(curInv, x0, r0, ci, cj, Aij, Aji)
    cij = (1-x0(ci)-r0(ci))*Aij;
    LUpdateij = curInv(:, ci);
    RUpdateij = curInv(cj, :);
    numeratorij = 1 + cij * curInv(cj,ci);
    curInv = curInv - cij* 1/numeratorij * (LUpdateij * RUpdateij);
    cji = (1-x0(cj)-r0(cj))*Aji;
    LUpdateji = curInv(:, cj);
    RUpdateji = curInv(ci, :);
    numeratorji = 1 + cji * curInv(ci,cj);
    udtInv = curInv - cji * 1/numeratorji * (LUpdateji * RUpdateji);
end