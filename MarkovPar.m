function [h] = MarkovPar(A,K,G,x)
%% System has n states, m inputs, l outputs
close all
l = size(G,1);
h = zeros(l,l,x+1);
title('n = 1')
pause(1)
for i = 1:x
    h(:,:,x) = G*(A-K*G)^(i-2)*K;
    imagesc(h(:,:,x));
    title(sprintf('n = %.0f',i))
    pause(1)
end
end