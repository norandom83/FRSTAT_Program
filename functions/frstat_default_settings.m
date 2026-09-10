function S = frstat_default_settings()
    S.bmgrid = 199;
    S.stepsize = 5;
    S.q_low = 0.05;
    S.q_high = 0.95;
    S.nRep = 1000;
    S.ciAlpha = 0.025;
    S.testReplications = 5000;
    S.testMcGrid = 499;
    S.loc_u_band = [0.25 0.75];
    S.loc_u_n = 201;
    S.loss_h = 0.50;
    S.loss_step = 0.05;
end
