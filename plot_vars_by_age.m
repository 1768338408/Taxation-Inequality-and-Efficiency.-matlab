function [fn] = plot_vars_by_age(fn, avgs_by_age, p)

n_age = p.n_age;


fn = fn + 1; figure(fn); clf;
plot([1:n_age], avgs_by_age.cons, 'b-', 'LineWidth', 2)
xlabel('age')
ylabel('consumption')
%print -dpng fig_C.png

fn = fn + 1; figure(fn); clf;
plot([1:n_age], avgs_by_age.K, 'b-', 'LineWidth', 2)
xlabel('age')
ylabel('capital')
%print -dpng fig_K.png

fn = fn + 1; figure(fn); clf;
plot([1:n_age], avgs_by_age.Kn, 'b-', 'LineWidth', 2)
xlabel('age')
ylabel('capital next-period')
%print -dpng fig_Kn.png

fn = fn + 1; figure(fn); clf;
plot([1:n_age], avgs_by_age.earn, 'b-', 'LineWidth', 2)
xlabel('age')
ylabel('earnings')
%print -dpng fig_earn.png
    
fn = fn + 1; figure(fn); clf;
plot([1:n_age], avgs_by_age.inc, 'b-', 'LineWidth', 2)
xlabel('age')
ylabel('income')
%print -dpng fig_inc.png

end