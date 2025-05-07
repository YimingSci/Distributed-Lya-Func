clear; clc;

% Initialize the network
Num_gen = 30;
K = 1;
Ks = 1.5;
Num_node = Num_gen;
mpc = WattsStrogatz(Num_node,1,0);
mpc.Edges.Weight = normrnd(-1, 0, size(mpc.Edges, 1), 1);
num_mpc = size(mpc.Edges,1);
J = zeros(Num_node);
for i_h = 1:num_mpc
    i_h1 = mpc.Edges{i_h,1}(1,1);
    i_h2 = mpc.Edges{i_h,1}(1,2);
    J(i_h1,i_h2) = mpc.Edges.Weight(i_h);
    J(i_h2,i_h1) = mpc.Edges.Weight(i_h);
end
load('randomEP.mat')
Matrix_sep = equilibria(:,1)';
Sto_acc = [];
sep_count = 0;

for i_sep = 1:1:size(Matrix_sep,1)
    i_sep = 1;
    phi_sep = Matrix_sep(i_sep,:);
    phi_star = phi_sep;
    J_matrix = zeros(Num_node, Num_node);
    for i = 1:Num_node
        for j = 1:Num_node
            if i == j
                sum_term = sum(J(i, [1:i-1, i+1:end]) .* cos(phi_star(i) - phi_star([1:i-1, i+1:end])));
                J_matrix(i,i) = -K * sum_term - 2 * Ks * cos(2 * phi_star(i));
            else
                J_matrix(i,j) = K * J(i,j) * cos(phi_star(i) - phi_star(j));
            end
        end
    end
    eig_out = max(real(eig(J_matrix)))


    if eig_out < 0
        Mat_Eta = [];
        Cell_V0 = [];
        Mat_lamda = [];
        Mat_w = [];
        Sto_B1_error_value = [];

        for i_wai = 1:1:10
            % Initialize the dynamical equations
            Zz = cell(2*Num_gen,1);
            for i = 1:1:2*Num_gen
                Zz{i,1} = pvar(['Z' num2str(i)]);
            end
            Z = [];
            for i = 1:1:2*Num_gen
                Z = [Z;Zz{i,1}];
            end
            for i = 1:1:Num_gen
                comp = 0;
                comp_2 = -2*Ks*(  Z(2*i-1)*cos(phi_sep(i)) + (1-Z(2*i))*sin(phi_sep(i)) ) * ( (1-Z(2*i))*cos(phi_sep(i)) - Z(2*i-1)*sin(phi_sep(i)) );
                for j = 1:1:Num_gen
                    if i == j
                    else
                        it1 = cos(phi_sep(i)-phi_sep(j));
                        it2 = sin(phi_sep(i)-phi_sep(j));
                        comp_1 = -K*J(i,j)*(it1*(Z(2*i-1)*(1-Z(2*j)) - Z(2*j-1)*(1-Z(2*i))) + it2*((1-Z(2*i))*(1-Z(2*j)) + Z(2*i-1)*Z(2*j-1)));
                        comp = comp + comp_1;
                    end
                end
                dtheta_dt = comp + comp_2;
                dZz_dt{2*i-1,1} = (1-Z(2*i))*dtheta_dt;
                dZz_dt{2*i,1} = Z(2*i-1)*dtheta_dt;
            end
            dZ_dt = [];
            for i = 1:1:2*Num_gen
                dZ_dt = [dZ_dt; dZz_dt{i,1}];
            end
            % Initialize the Lyapunov function
            if i_wai == 1
                V0 = 0;
                for i = 1:1:Num_gen
                    V0 = V0 +  ((pi/2 - ( 1 - (Z(2*i-1)^2 + Z(2*i)^2)/2 )) + phi_sep(i))^2;
                end
                Eta = 100;
                V0 = V0/Eta;
                Vdot_z = jacobian(V0, Z);
                Vdot0 = Vdot_z * dZ_dt;
                Eta = 1;
            end

            pvar gamma
            G = [];
            for i = 1:1:Num_gen
                G = [G; Z(2*i-1)^2 + Z(2*i)^2 - 2*Z(2*i)];
            end

            eps = 0.01;
            kst = 0; 
            lf = Num_gen;
            
            Cell_Zdec1 = cell(lf,1);
            Cell_Zrest1 = cell(lf,1);
            Cell_dZdt_dec1 = cell(lf,1);
            Cell_p0 = cell(lf,1);
            Cell_dp0_dt = cell(lf,1);
            Cell_indices_1 = cell(lf,1);
            Cell_Zdec2 = cell(lf,1);
            Cell_Zrest2 = cell(lf,1);
            Cell_dZdt_dec2 = cell(lf,1);
            Cell_v0 = cell(lf,1);
            Cell_vdot0 = cell(lf,1);
            Cell_indices_2 = cell(lf,1);
            Matrix_Eta_2 = ones(lf,1);
            % Initialize V_{k,\ell}
            for i_ci = 1:1:lf
                Z_dec1 = [];
                indices_1 = mod([i_ci,i_ci+1]- 1, Num_gen) + 1;
                sorted_cdc_indices = sort(indices_1);
                for i_cdc = 1:1:size(sorted_cdc_indices,2)
                    Z_dec1 = [Z_dec1; Zz{2*(sorted_cdc_indices(i_cdc))-1,1};Zz{2*(sorted_cdc_indices(i_cdc)),1}];
                end
                Z_rest1 = [];
                for i_cdc = 1:Num_gen
                    if ~ismember(i_cdc, indices_1)
                        Z_rest1 = [Z_rest1; Z(2*i_cdc-1); Z(2*i_cdc)];
                    end
                end

                dZ_dt_dec_pre1 = subs(dZ_dt,Z_rest1,[zeros(size(Z_rest1,1),1)]);
                dZ_dt_dec1 = [];
                for i_cdc = 1:1:size(indices_1,2)
                    dZ_dt_dec1 = [dZ_dt_dec1;dZ_dt_dec_pre1(2*sorted_cdc_indices(i_cdc)-1,1);dZ_dt_dec_pre1(2*sorted_cdc_indices(i_cdc),1)];
                end

                Z_rest1_0 = zeros(size(Z_rest1,1),1);
                p0 = subs(V0, Z_rest1, Z_rest1_0);
                pdot_z = jacobian(p0, Z_dec1);
                dp0_dt = pdot_z * dZ_dt_dec1;
                Z_dec2 = [];

                if i_ci == Num_gen-1
                    indices_2 = [Num_gen-1, Num_gen, 1];
                elseif i_ci == Num_gen
                    indices_2 = [Num_gen, 1, 2];
                else
                    indices_2 = [i_ci, i_ci+1, i_ci+2];
                end
                sorted_cdc_indices = sort(indices_2);
                for i_cdc = 1:1:size(sorted_cdc_indices,2)
                    Z_dec2 = [Z_dec2; Zz{2*(sorted_cdc_indices(i_cdc))-1,1};Zz{2*(sorted_cdc_indices(i_cdc)),1}];
                end
                Z_rest2 = [];
                for i_cdc = 1:Num_gen
                    if ~ismember(i_cdc, indices_2)
                        Z_rest2 = [Z_rest2; Z(2*i_cdc-1); Z(2*i_cdc)];
                    end
                end
                dZ_dt_dec_pre2 = subs(dZ_dt,Z_rest2,[zeros(size(Z_rest2,1),1)]);
                dZ_dt_dec2 = [];
                for i_cdc = 1:1:size(indices_2,2)
                    dZ_dt_dec2 = [dZ_dt_dec2;dZ_dt_dec_pre2(2*sorted_cdc_indices(i_cdc)-1,1);dZ_dt_dec_pre2(2*sorted_cdc_indices(i_cdc),1)];
                end
                Z_rest2_0 = zeros(size(Z_rest2, 1),1);
                v0 = subs(V0, Z_rest2, Z_rest2_0);
                vdot_z = jacobian(v0, Z_dec2);
                vdot0 = vdot_z * dZ_dt_dec2;
                Eta_2 = 1;
                Cell_Zdec1{i_ci,1} = Z_dec1;
                Cell_Zrest1{i_ci,1} = Z_rest1;
                Cell_dZdt_dec1{i_ci,1} = dZ_dt_dec1;
                Cell_p0{i_ci,1} = p0;
                Cell_dp0_dt{i_ci,1} = dp0_dt;
                Cell_indices_1{i_ci,1} = indices_1;
                Cell_Zdec2{i_ci,1} = Z_dec2;
                Cell_Zrest2{i_ci,1} = Z_rest2;
                Cell_dZdt_dec2{i_ci,1} = dZ_dt_dec2;
                Cell_v0{i_ci,1} = v0;
                Cell_vdot0{i_ci,1} = vdot0;
                Cell_indices_2{i_ci,1} = indices_2;
            end
            Cell_alpha_ca = cell(lf,size(Cell_indices_2{i_ci,1},2) - 1);
            Cell_beta_ca = cell(lf,size(Cell_indices_2{i_ci,1},2) - 1);
            Cell_tau_ca = cell(lf,size(Cell_indices_2{i_ci,1},2) - 1);
            gamma_mat_var = 0.5*ones(lf,1);


            for i_ci = 1:1:lf
                Store_A1_gamma = [];
                Sto_A1_error_value = [];
                Store_p = [];
                Sto_error_value = [];
                for i_A1 = 1:1:10
                    % Algorithm 2: SOS-based estimation of the optimal ROA boundary
                    s1_deg = monomials(Cell_Zdec1{i_ci,1}, 1);
                    Lst = 0.001*(Cell_Zdec2{i_ci,1}).' * Cell_Zdec2{i_ci,1};
                    iter = 0;
                    feasibility = false;
                    while_count = 0;
                    gamma_low = 0; gamma_high = 10;
                    while gamma_high - gamma_low > eps || ~feasibility
                        gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
                        disp("----------");
                        disp("check gamma: ");
                        disp(gamma_var);
                        gamma_mat_var(i_ci,1) = gamma_var;
                        prog_A1 = sosprogram(Cell_Zdec2{i_ci,1});
                        s11_deg = monomials(Cell_Zdec1{i_ci,1}, 0);
                        s2_deg = monomials(Cell_Zdec2{i_ci,1}, 1);
                        lamda_deg = monomials(Cell_Zdec1{i_ci,1}, 0);
                        L1 = 1e-6 * (Cell_Zdec1{i_ci,1}).' * (Cell_Zdec1{i_ci,1});
                        L2 = L1; L3 = L1;
                        [prog_A1, s1] = sossosvar(prog_A1, s1_deg);
                        [prog_A1, s11] = sossosvar(prog_A1, s11_deg);
                        [prog_A1, s2] = sossosvar(prog_A1, s2_deg);
                        [prog_A1, s4] = sossosvar(prog_A1, s2_deg);
                        [prog_A1, s41] = sossosvar(prog_A1, s11_deg);
                        [prog_A1, lamda11] = sossosvar(prog_A1, lamda_deg);
                        [prog_A1, lamda12] = sossosvar(prog_A1, lamda_deg);
                        [prog_A1, lamda21] = sossosvar(prog_A1, lamda_deg);
                        [prog_A1, lamda22] = sossosvar(prog_A1, lamda_deg);
                        prog_A1 = sosineq(prog_A1, s1*(Cell_p0{i_ci,1}-gamma_var+kst*Lst) - s11*(Cell_v0{i_ci,1}-Matrix_Eta_2(i_ci,1)) - lamda21*G(Cell_indices_1{i_ci,1}(1)) - lamda22*G(Cell_indices_1{i_ci,1}(2)) );
                        prog_A1 = sosineq(prog_A1, Cell_p0{i_ci,1} + kst*Lst - L1 - lamda11*G(Cell_indices_1{i_ci,1}(1)) - lamda12*G(Cell_indices_1{i_ci,1}(2)) );
                        prog_A1 = sosineq(prog_A1, (Cell_v0{i_ci,1} - Matrix_Eta_2(i_ci,1)) * s4  - s41*(L2 + Cell_vdot0{i_ci,1} ) );
                        tau_deg = monomials(Cell_Zdec1{i_ci,1}, 1);
                        alpha_deg = monomials(Cell_Zdec1{i_ci,1}, 0);
                        beta_deg = monomials(Cell_Zdec1{i_ci,1}, 0);
                        for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                            [prog_A1, Cell_alpha_ca{i_ci,ioo_nei}] = sossosvar(prog_A1, alpha_deg);
                            [prog_A1, Cell_beta_ca{i_ci,ioo_nei}] = sossosvar(prog_A1, beta_deg);
                            [prog_A1, Cell_tau_ca{i_ci,ioo_nei}] = sossosvar(prog_A1, tau_deg);
                        end
                        Sto_ca = 0;
                        for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                            indice_vir = Cell_indices_2{i_ci,1}(1,ioo_nei);
                            Sto_ca = Sto_ca + Cell_p0{indice_vir,1}*(Cell_alpha_ca{i_ci,ioo_nei}...
                                - Cell_beta_ca{i_ci,ioo_nei})...
                                + Cell_beta_ca{i_ci,ioo_nei}*Cell_v0{i_ci,1}...
                                + Cell_tau_ca{i_ci,ioo_nei}*(Cell_p0{indice_vir,1}-gamma_mat_var(indice_vir,1));
                        end
                        prog_A1 = sosineq(prog_A1,  -Cell_dp0_dt{i_ci,1} - L3 + Sto_ca );

                        ss_AB = 0;
                        ss_B = 0;
                        for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                            if ioo_nei == 1
                            else
                                prog_A1 = sosineq(prog_A1,  (Cell_alpha_ca{i_ci,ioo_nei} - Cell_beta_ca{i_ci,ioo_nei}) - (1e-2) );
                            end
                            ss_AB = ss_AB + (Cell_alpha_ca{i_ci,ioo_nei} - Cell_beta_ca{i_ci,ioo_nei});
                            ss_B = ss_B + (Cell_beta_ca{i_ci,ioo_nei});
                        end
                        prog_A1 = sosineq(prog_A1,  -(ss_AB) - (1e-2) );
                        prog_A1 = sosineq(prog_A1,  (ss_B) - (1e-2) );


                        prog_A1 = sosineq(prog_A1, (Cell_p0{i_ci,1} - gamma_var) * s2  - (L2 + Cell_dp0_dt{i_ci,1} ));
                        prog_A1 = sossolve(prog_A1);
                        if prog_A1.solinfo.info.pinf || prog_A1.solinfo.info.dinf
                            disp("Infeasible");
                            gamma_high = gamma_var;
                            feasibility = false;
                        elseif prog_A1.solinfo.info.numerr == 1
                            disp("Numerical_error");
                            gamma_high = gamma_high - 0.5;
                            feasibility = false;
                        else
                            disp("Feasible");
                            gamma_low = gamma_var;
                            feasibility = true;
                        end
                        if while_count > 60
                            disp("Infeasible in B1-1");
                            break
                        end
                        while_count = while_count +1;
                    end
                    gamma_in = gamma_var;
                    s1_in = sosgetsol(prog_A1,s1);
                    s2_in = sosgetsol(prog_A1,s2);
                    s4_in = sosgetsol(prog_A1,s4);
                    s11_in = sosgetsol(prog_A1,s11);
                    s41_in = sosgetsol(prog_A1,s41);
                    lamda11_in = sosgetsol(prog_A1,lamda11);
                    lamda12_in = sosgetsol(prog_A1,lamda12);
                    lamda21_in = sosgetsol(prog_A1,lamda21);
                    lamda22_in = sosgetsol(prog_A1,lamda22);
                    for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                        Cell_alpha_in{i_ci,ioo_nei} = double(sosgetsol(prog_A1,Cell_alpha_ca{i_ci,ioo_nei}));
                        Cell_beta_in{i_ci,ioo_nei} = double(sosgetsol(prog_A1,Cell_beta_ca{i_ci,ioo_nei}));
                        Cell_tau_in{i_ci,ioo_nei} = sosgetsol(prog_A1,Cell_tau_ca{i_ci,ioo_nei});
                    end
                    % Algorithm 3: Construction of the local Lyapunov function
                    prog_A4 = sosprogram(Cell_Zdec2{i_ci,1});
                    p_deg = monomials(Cell_Zdec1{i_ci,1}, 0:1);
                    [prog_A4, p0_var] = sossosvar(prog_A4, p_deg);
                    pdot_z_var = jacobian(p0_var, Cell_Zdec1{i_ci,1});
                    pdot0_var = pdot_z_var * Cell_dZdt_dec1{i_ci,1};
                    prog_A4 = sosineq(prog_A4, s1_in*(p0_var-gamma_in+kst*Lst) - ...
                        s11_in*(Cell_v0{i_ci,1}-Matrix_Eta_2(i_ci,1)) - lamda21_in...
                        *G(Cell_indices_1{i_ci,1}(1))...
                        - lamda22_in*G(Cell_indices_1{i_ci,1}(2)) );
                    prog_A4 = sosineq(prog_A4, p0_var + kst*Lst - L1 - lamda11_in*...
                        G(Cell_indices_1{i_ci,1}(1))...
                        - lamda12_in*G(Cell_indices_1{i_ci,1}(2)) );
                    prog_A4 = sosineq(prog_A4, (Cell_v0{i_ci,1} - Matrix_Eta_2(i_ci,1)) * s4_in...
                        - s41_in*(L2 + Cell_vdot0{i_ci,1} ) );
                    Sto_in = 0;
                    for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                        indice_vir = Cell_indices_2{i_ci,1}(1,ioo_nei);
                        if ioo_nei == 1
                            Sto_in = Sto_in + p0_var*(Cell_alpha_in{i_ci,ioo_nei}...
                                - Cell_beta_in{i_ci,ioo_nei})...
                                + Cell_beta_in{i_ci,ioo_nei}*Cell_v0{i_ci,1}...
                                + Cell_tau_in{i_ci,ioo_nei}*(p0_var-gamma_mat_var(indice_vir,1));
                        else
                            Sto_in = Sto_in + Cell_p0{indice_vir,1}*(Cell_alpha_in{i_ci,ioo_nei}...
                                - Cell_beta_in{i_ci,ioo_nei})...
                                + Cell_beta_in{i_ci,ioo_nei}*Cell_v0{i_ci,1}...
                                + Cell_tau_in{i_ci,ioo_nei}*(Cell_p0{indice_vir,1}-...
                                gamma_mat_var(indice_vir,1));
                        end
                    end
                    prog_A4 = sosineq(prog_A4,  -pdot0_var - L3 + Sto_in );
                    ss_AB = 0;
                    ss_B = 0;
                    for ioo_nei = 1:1:size(Cell_indices_2{i_ci,1},2) - 1
                        if ioo_nei == 1
                        else
                            prog_A4 = sosineq(prog_A4,  (Cell_alpha_in{i_ci,ioo_nei} - Cell_beta_in{i_ci,ioo_nei}) - (1e-2) );
                        end
                        ss_AB = ss_AB + (Cell_alpha_in{i_ci, ioo_nei} - Cell_beta_in{i_ci,ioo_nei});
                        ss_B = ss_B + (Cell_beta_in{i_ci, ioo_nei});
                    end
                    prog_A4 = sosineq(prog_A4,  -(ss_AB) - (1e-2) );
                    prog_A4 = sosineq(prog_A4,  (ss_B) - (1e-2) );
                    prog_A4 = sosineq(prog_A4, (p0_var - gamma_in) * s2_in  - (L2 + pdot0_var ));
                    prog_A4 = sossolve(prog_A4);
                    p0_in = sosgetsol(prog_A4, p0_var);
                    Cell_p0{i_ci,1} = p0_in/gamma_mat_var(i_ci,1);
                    Cell_dp0_dt{i_ci,1} = jacobian(Cell_p0{i_ci,1}, Cell_Zdec1{i_ci,1})...
                        * Cell_dZdt_dec1{i_ci,1};
                    Store_A1_gamma = [Store_A1_gamma; gamma_var];
                    Store_p = [Store_p;p0_in/gamma_mat_var(i_ci,1)];


                    if i_A1 > 1
                        p_error = Store_p(end) - Store_p(end-1);
                        p_error_value = double(max(p_error.coefficient+0));
                        Sto_error_value = [Sto_error_value;p_error_value]
                        if double(max(p_error.coefficient+0)) < 0.1
                            break;
                        end
                    end
                end
            end

            % Construct the global Lyapunov function
            A = [];
            B = [];
            max_col = lf;  
            for i = 1:lf
                for j = 1:lf
                    if j == i
                        A(i,j) = double(Cell_alpha_in{i,1});
                        B(i,j) = double(Cell_beta_in{i,1});
                    elseif j == i+1
                        A(i,j) = double(Cell_alpha_in{i,2});
                        B(i,j) = double(Cell_beta_in{i,2});
                    else
                        A(i,j) = 0;
                        B(i,j) = 0;
                    end

                end
            end

            Out1 = A-B;
            Outt = sum(Out1,2);

            StoAB = Out1;  
            StoB = B;    

            lf = size(StoAB, 1);

            nvars = lf + 1;  
            w0 = [ones(lf,1)/lf; 0.1]; 

            obj_fun = @(x) -x(end);  

            Aeq = [ones(1, lf), 0]; 
            beq = 1;  

            A = [];  
            b = [];  

            nonlcon = @(x) myNonlcon2(x, StoAB, StoB);  

            epsilon = 0;  
            lb = [epsilon * ones(lf,1); 0.1]; 
            ub = [(1 - epsilon) * ones(lf,1); 1]; 
            options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp'); 
            [x_opt, fval, exitflag] = fmincon(obj_fun, w0, A, b, Aeq, beq, lb, ub, nonlcon, options);  
            w_matrix = x_opt(1:lf); 
            Vg = 0;
            for iop = 1:1:lf
                Vg = Vg + Cell_p0{iop,1}*w_matrix(iop,1);
            end
            Eta = (gamma_mat_var')*w_matrix;
            Mat_Eta = [Mat_Eta;Eta];

            V0 = Vg/Eta;
            Cell_V0 = [Cell_V0; V0];

            Mat_lamda = [Mat_lamda; x_opt(end)];
            Mat_w =[Mat_w, w_matrix];

            if i_wai > 1
                V_error = Cell_V0(end) - Cell_V0(end-1);
                V_error_value = double(max(V_error.coefficient+0));
                Sto_V_error_value = [Sto_V_error_value;V_error_value]
                if double(max(V_error.coefficient+0)) < eps
                    break;
                end
            end
        end
    end
end
