function f_ref0 = build_reference_density_seed(y_grid, yvar0)
% Build reference density seed ONLY from yvar0 (temperature PDF).

    y_grid = y_grid(:);

    if ismatrix(yvar0) && size(yvar0,2) > 1
        f_ref0 = mean(yvar0,2,'omitnan');
    else
        f_ref0 = yvar0(:);
    end

    f_ref0 = max(f_ref0,0);

    area0 = trapz(y_grid, f_ref0);
    if area0 <= 0
        error('Reference density seed has non-positive integral. Check yvar0 and y_grid.');
    end

    f_ref0 = f_ref0 / area0;
end