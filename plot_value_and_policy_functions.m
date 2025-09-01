function [fn] = plot_value_and_policy_functions(fn, a_grid, v, g, p)

n_age = p.n_age;

age_vec = [1, 30, n_age-1, n_age];
for i = 1:length(age_vec)
    age = age_vec(i);

    fn = fn + 1; figure(fn); clf;
    plot(a_grid, v(:,:,age), 'b-', 'LineWidth', 2)
    grid on
    title(sprintf('value functions, age = %d', age))
    xlabel('assets')

    fn = fn + 1; figure(fn); clf;
    hold on
    plot(a_grid, g(:,:,age), 'b-', 'LineWidth', 2)
    plot(a_grid, a_grid, 'k', 'LineWidth', 2)
    hold off
    grid on

end

end