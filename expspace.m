function [x_vec, s_x] = expspace(x_min, x_max, x_scale, n_x)

    if (n_x > 1)
        idx   = linspace(1, n_x, n_x)';
        s_x   = ( (idx-1)./(n_x-1) ).^x_scale;
        x_vec = x_min + (x_max-x_min).*s_x;
    else
        s_x = 1;
        x_vec = 0.5*(x_min + x_max);
    end

end