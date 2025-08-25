function generate_experiment_matrices(m, n, r_values, d_values, n_mtx, output_dir)
    % Creates output directory if it does not exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % seed = randi([1, 1000000]);
    seed = 2;

    rng(seed);

    % Iterates over all values of r and d
    for i = 1:length(r_values)
        for j = 1:length(d_values)
            r = r_values(i);
            d = d_values(j);
            
            % Generates n_mtx matrices for each combination of r and r
            for k = 1:n_mtx
                % Generates matrix using function generate_experiment_matrix
                A = gen_single_matrix(m, n, r, d);
                
                % Creates matrix file name
                d_decimal = round(d * 100);
                filename = sprintf('A_m%d_n%d_r%d_d%d_idx%d.mat', m, n, r, d_decimal, k);
                filepath = fullfile(output_dir, filename);
                
                % Saves matrix as a file .mat
                save(filepath, 'A');
                fprintf('Salved: %s', filepath);
            end
        end
    end
end