function L = truncLoss(L0,tao)

L0(abs(L0)>tao) = tao;

L = L0;