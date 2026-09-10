function W = brownian_motion(t0, T, N, x0)
    dt = (T - t0) / N;
    dW = sqrt(dt) * randn(1, N);
    W  = [x0, x0 + cumsum(dW)];
end