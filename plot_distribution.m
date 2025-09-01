function [fn] = plot_distribution(fn, a_grid2, Phi, Phi_all, p)

n_age = p.n_age;

Phi_all_sum = sum(Phi_all,2);

fn = fn + 1; figure(fn); clf;
plot(a_grid2, Phi_all_sum, 'b-', 'LineWidth', 2);
title('overall cdf')
xlabel('assets')
ylabel('fraction of agents')


age_vec = [1, 30, n_age-1, n_age];
for i = 1:length(age_vec)
    age = age_vec(i);

    Phi_age = sum(Phi(:,:,age),2);
    fn = fn + 1; figure(fn); clf;
    plot(a_grid2, Phi_age, 'b-', 'LineWidth', 2)
    grid on
    title(sprintf('cdf, age = %d', age))
    xlabel('assets')
    ylabel('fraction of agents')


end

end