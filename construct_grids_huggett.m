function [ a_grid, z_grid, a_grid2, Pz, Pz_cdf, Pz1, Pz1_cdf ] = construct_grids_huggett( p )

% get number of grid dpoints
n_a   = p.n_a;
n_z   = p.n_z;
n_a2  = p.n_a2;

% grid points for assets
a_grid = expspace(p.a_min, p.a_max, p.a_scale, n_a);

% grid pints for assets (more fine)
a_grid2 = expspace(p.a_min, p.a_max, p.a_scale, n_a2);

% grid points for endowment shock, z
extra_grid_point = p.z_grid_opt == 2;
z_min = p.z_min_w*p.sd_y1;
z_max1 = p.z_max_w1*p.sd_y1;
z_max2 = p.z_max_w2*p.sd_y1;

if (extra_grid_point)
    z_grid1 = linspace(z_min, z_max1, n_z-1)';
    z_grid = [z_grid1; z_max2];
else
    z_grid = linspace(z_min, z_max1, n_z);
end

% construct midpoint of grid points for earnings
midpoint = 0.5*(z_grid(2:end) + z_grid(1:end-1));

% construct tauchen approximation for earnings process
Pz = zeros(n_z,n_z);
for i = 1:n_z
    Pz(i,1)   = normcdf( (midpoint(1) - p.gamma*z_grid(i) ) /  p.sd_e);
    Pz(i,n_z) = normcdf( (midpoint(n_z-1) - p.gamma*z_grid(i) ) /  p.sd_e , 'upper');

    % compute probabilities for remaining grid points
    for j = 2:n_z-1
        Pz(i,j) = normcdf( ( midpoint(j) - p.gamma*z_grid(i) ) / p.sd_e ) ...
            - normcdf( ( midpoint(j-1) - p.gamma*z_grid(i) ) / p.sd_e );
    end
end


% construct cdf for probability transition matrix
Pz_cdf = cumsum( Pz, 2 );

% construct tauchen approximation for initial earnings
Pz1 = zeros(p.n_z, 1);

% using false input, alnorm returns normal cdf
Pz1(1) = normcdf( midpoint(1) / p.sd_y1 );

% using false input, alnorm returns 1 - cdf
Pz1(n_z) = normcdf( midpoint(n_z-1) / p.sd_y1 , 'upper' );

% compute probabilities for remaining grid points
for j = 2:p.n_z-1
    Pz1(j) = normcdf( midpoint(j) / p.sd_y1) ...
        - normcdf( midpoint(j-1) / p.sd_y1 );
end

% construct cdf for probability transition matrix
Pz1_cdf = cumsum( Pz1 );
    
end

