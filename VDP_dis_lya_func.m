clear; clc;

% Initialize the network
Num_gen = 5;
Sto_zong = [];
Alpha_mean = 1.5;
Alpha_var = Alpha_mean * 0.2;

for i_alpha = 1:1:100
    Alpha = normrnd(Alpha_mean, Alpha_var, [Num_gen, 1]);
    if min(Alpha) > 0
        break
    end
end

i_iter1_max = 100;
for i_iter1 = 1:1:i_iter1_max

    K_mean = 0.02;
    K_var = 0.0001;
    mpc = WattsStrogatz(Num_gen,2,0);
    mpc.Edges.Weight = normrnd(K_mean, K_var, size(mpc.Edges, 1), 1);

    num_mpc = size(mpc.Edges,1);
    K_init = zeros(Num_gen);
    K = zeros(Num_gen);
    for i_h = 1:num_mpc
        i_h1 = mpc.Edges{i_h,1}(1,1);
        i_h2 = mpc.Edges{i_h,1}(1,2);
        K_init(i_h1,i_h2) = 1;
        K_init(i_h2,i_h1) = 1;
        K(i_h1,i_h2) = mpc.Edges.Weight(i_h);
        K(i_h2,i_h1) = mpc.Edges.Weight(i_h);
    end

    % Initialize the dynamical equations
    y = cell(Num_gen,1);
    x = cell(Num_gen,1);
    for i = 1:1:Num_gen
        y{i,1} = pvar(['y' num2str(i)]);
        x{i,1} = pvar(['x' num2str(i)]);
    end
    Z = [];
    for i = 1:1:Num_gen
        Z = [Z;y{i,1};x{i,1}];
    end

    for i = 1:1:Num_gen
        dy_dt{i,1} = - Z(2*i);
        Coup_mat1 = 0;
        for j = 1:1:Num_gen
            if i ~= j
                Coup_mat1 = Coup_mat1 + K(i,j)*(Z(2*j-1)-Z(2*i-1));
            end
        end
        dx_dt{i,1} = Z(2*i-1) + Alpha(i)*Z(2*i)*(-1+Z(2*i-1)^2) + Coup_mat1;
    end
    dZ_dt = [];
    for i = 1:1:Num_gen
        dZ_dt = [dZ_dt; dy_dt{i,1}; dx_dt{i,1}];
    end

    % Initialize the Lyapunov function
    jac = jacobian(dZ_dt, Z);
    A = double(subs(jac,Z,zeros(2*Num_gen,1)));
    Q = eye(2*Num_gen);
    P = lyap(A', Q);
    V0 = Z' * P * Z;
    Vdot_z = jacobian(V0, Z);
    Vdot0 = Vdot_z * dZ_dt;

    Cell_V0 = [];
    Mat_Eta = [];
    Mat_lamda = [];

    Cell_V0 = [Cell_V0; V0];
    Eta = 3;

    Mat_Eta = [Mat_Eta; Eta];
    Mat_w = [];
    Sto_V_error_value = [];
    V0 = V0/Eta;

    for i_wai = 1:1:10
        Vdot_z = jacobian(V0, Z);
        Vdot0 = Vdot_z * dZ_dt;
        Eta = 1;
        V_jac = V0;
        Eta_jac = 1;

        pvar gamma
        k = 2; 
        lf = factorial(Num_gen) / (factorial(k) * factorial(Num_gen-k)) +Num_gen;

        Cell_Zdec = cell(lf,1);
        Cell_Zdec_rest = cell(lf,1);
        Cell_Vdec = cell(lf,1);
        Cell_dVdec_dt = cell(lf,1);
        Cell_dZdec_dt = cell(lf,1);

        % Initialize V_{k,\ell}
        i_cout = 1;
        for i = 1:1:Num_gen
            for j = i:1:Num_gen
                if i == j
                    Cell_Zdec{i_cout,1} = [y{i,1};x{i,1}];
                else
                    Cell_Zdec{i_cout,1} = [y{i,1};x{i,1};y{j,1};x{j,1}];
                end
                Z_rest = [];
                for it = 1:1:Num_gen
                    if it ~= i
                        if it ~= j
                            Z_rest = [Z_rest;y{it,1};x{it,1}];
                        end
                    end
                end
                Cell_Zdec_rest{i_cout,1} = Z_rest;
                Cell_Vdec{i_cout,1} = subs(V0,Z_rest,[zeros(size(Z_rest,1),1)]);
 
                Cell_dZdec_dt{i_cout,1} = subs(dZ_dt,Z_rest,[zeros(size(Z_rest,1),1)]);
                Cell_dVdec_dt{i_cout,1} = jacobian(Cell_Vdec{i_cout,1}, Z)*Cell_dZdec_dt{i_cout,1};

                i_cout = i_cout + 1;
            end
        end
        eps = 1e-3;
        s1_deg = monomials(Z, 0);
        Lst = 0.5*(Z).' * Z;
        kst = 0.01;
        gamma_mat_var = 0.5*ones(lf,1);
        i_iter2_max = 5;
        Cell_alpha_ca = cell(lf,lf);
        Cell_beta_ca = cell(lf,lf);
        Cell_tau_ca = cell(lf,lf);
        Cell_alpha_in = cell(lf,lf);
        Cell_beta_in = cell(lf,lf);
        Cell_tau_in = cell(lf,lf);
        for ioo = 1:1:lf
            Store_A1_gamma = [];
            Sto_A1_error_value = [];
            Store_p = [];
            Sto_error_value = [];
            for i_iter2 = 1:1:i_iter2_max
                % Algorithm 2: SOS-based estimation of the optimal ROA boundary
                iter = 0;
                feasibility = false;
                while_count = 0;
                gamma_low = 0; gamma_high = 5;
                while gamma_high - gamma_low > eps || ~feasibility
                    gamma_var = (gamma_high - gamma_low) / 2 + gamma_low;
                    disp("----------");
                    disp("check gamma: ");
                    disp(gamma_var);
                    gamma_mat_var(ioo,1) = gamma_var;
                    prog_A1 = sosprogram(Z);
                    s11_deg = monomials(Cell_Zdec{ioo,1}, 0);
                    s2_deg = monomials(Cell_Zdec{ioo,1}, 1);
                    s3_deg = monomials(Z, 1);
                    s31_deg = monomials(Z, 0);
                    tau_deg = monomials(Cell_Zdec{ioo,1}, 1);
                    alpha_deg = monomials(Cell_Zdec{ioo,1}, 0);
                    beta_deg = monomials(Cell_Zdec{ioo,1}, 0);
                    L1 = eps * (Cell_Zdec{ioo,1}).' * (Cell_Zdec{ioo,1});
                    L2 = L1;
                    L3 = L1;                 
                    [prog_A1, s1] = sossosvar(prog_A1, s1_deg);
                    [prog_A1, s11] = sossosvar(prog_A1, s11_deg);
                    [prog_A1, s2] = sossosvar(prog_A1, s2_deg);
                    [prog_A1, s3] = sossosvar(prog_A1, s3_deg);
                    [prog_A1, s31] = sossosvar(prog_A1, s31_deg);
                    for ioo_nei = 1:1:lf
                        [prog_A1, Cell_alpha_ca{ioo,ioo_nei}] = sossosvar(prog_A1, alpha_deg);
                        [prog_A1, Cell_beta_ca{ioo,ioo_nei}] = sossosvar(prog_A1, beta_deg);
                        [prog_A1, Cell_tau_ca{ioo,ioo_nei}] = sossosvar(prog_A1, tau_deg);
                    end
                    prog_A1 = sosineq(prog_A1, Cell_Vdec{ioo,1} + kst*Lst - L1 );
                    prog_A1 = sosineq(prog_A1, s1*(Cell_Vdec{ioo,1}-gamma_var+kst*Lst) - s11*(V0-Eta) );
                    prog_A1 = sosineq(prog_A1, (Cell_Vdec{ioo,1} - gamma_var) * s2  - (L2 + Cell_dVdec_dt{ioo,1} ));
                    Sto_ca = 0;
                    for ioo_nei = 1:1:lf
                        Sto_ca = Sto_ca + Cell_Vdec{ioo_nei,1}*(Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei})...
                            + Cell_beta_ca{ioo,ioo_nei}*V0 + Cell_tau_ca{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}-gamma_mat_var(ioo_nei,1));
                    end
                    prog_A1 = sosineq(prog_A1,  -Cell_dVdec_dt{ioo,1} - L3 + Sto_ca );
                    prog_A1 = sosineq(prog_A1, (V0 - Eta) * s3  - s31*(L2 + Vdot0 ));  
                    ss_AB = 0;
                    ss_B = 0;
                    for ioo_nei = 1:1:lf
                        if ioo == ioo_nei                          
                        else 
                            prog_A1 = sosineq(prog_A1,  (Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei}) - (eps) );
                        end
                        ss_AB = ss_AB + (Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei});
                        ss_B = ss_B + (Cell_beta_ca{ioo,ioo_nei});
                    end
                    prog_A1 = sosineq(prog_A1,  -(ss_AB) - (eps) );
                    prog_A1 = sosineq(prog_A1,  (ss_B) - (eps) );                  
                    prog_A1 = sossolve(prog_A1);
                    if prog_A1.solinfo.info.pinf || prog_A1.solinfo.info.dinf
                        disp("Infeasible");
                        gamma_high = gamma_var;
                        feasibility = false;
                    elseif prog_A1.solinfo.info.numerr == 1
                        disp("Numerical_error");
                        gamma_high = gamma_high - eps;
                        feasibility = false;
                    else
                        disp("Feasible");
                        gamma_low = gamma_var;
                        feasibility = true;
                    end
                    if while_count > 100
                        disp("Infeasible in B1-1");
                        break
                    end
                    while_count = while_count +1;
                end
                gamma_var
                s1_in = sosgetsol(prog_A1,s1);
                s2_in = sosgetsol(prog_A1,s2);
                s11_in = sosgetsol(prog_A1,s11);
                s3_in = sosgetsol(prog_A1,s3);
                s31_in = sosgetsol(prog_A1,s31);
                for ioo_nei = 1:1:lf
                    Cell_alpha_in{ioo,ioo_nei} = double(sosgetsol(prog_A1,Cell_alpha_ca{ioo,ioo_nei}));
                    Cell_beta_in{ioo,ioo_nei} = double(sosgetsol(prog_A1,Cell_beta_ca{ioo,ioo_nei}));
                    Cell_tau_in{ioo,ioo_nei} = sosgetsol(prog_A1,Cell_tau_ca{ioo,ioo_nei});
                end
                prog_A2 = sosprogram(Z,gamma);
                prog_A2 = sosineq(prog_A2, Cell_Vdec{ioo,1} + kst*Lst - L1 );
                prog_A2 = sosineq(prog_A2, s1_in*(Cell_Vdec{ioo,1}-gamma+kst*Lst) - s11_in*(V0-Eta) );
                prog_A2 = sosineq(prog_A2, (Cell_Vdec{ioo,1} - gamma) * s2_in  - (L2 + Cell_dVdec_dt{ioo,1} ));
                Sto_in = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                        Sto_in = Sto_in + Cell_Vdec{ioo_nei,1}*(Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei})...
                            + Cell_beta_in{ioo,ioo_nei}*V0 + Cell_tau_in{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}- gamma);
                    else
                        Sto_in = Sto_in + Cell_Vdec{ioo_nei,1}*(Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei})...
                            + Cell_beta_in{ioo,ioo_nei}*V0 + Cell_tau_in{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}-gamma_mat_var(ioo_nei,1));
                    end
                end
                prog_A2 = sosineq(prog_A2,  -Cell_dVdec_dt{ioo,1} - L3 + Sto_in );
                prog_A2 = sosineq(prog_A2, (V0 - Eta) * s3_in  - s31_in*(L2 + Vdot0 ));
                ss_AB = 0;
                ss_B = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                    else
                        prog_A2 = sosineq(prog_A2,  (Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei}) - (eps) );
                    end
                    ss_AB = ss_AB + (Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei});
                    ss_B = ss_B + (Cell_beta_in{ioo,ioo_nei});
                end
                prog_A2 = sosineq(prog_A2,  -(ss_AB) - (eps) );
                prog_A2 = sosineq(prog_A2,  (ss_B) - (eps) );
                prog_A2 = sossetobj(prog_A2, -gamma);
                prog_A2 = sossolve(prog_A2);
                gamma_in = sosgetsol(prog_A2, gamma);
                gamma_mat_var(ioo,1) = gamma_in;

                % Algorithm 3: Construction of the local Lyapunov function
                pvar kst_var
                prog_A3 = sosprogram(Z,kst_var);
                [prog_A3, s1] = sossosvar(prog_A3, s1_deg);
                [prog_A3, s11] = sossosvar(prog_A3, s11_deg);
                [prog_A3, s2] = sossosvar(prog_A3, s2_deg);
                for ioo_nei = 1:1:lf
                    [prog_A3, Cell_alpha_ca{ioo,ioo_nei}] = sossosvar(prog_A3, alpha_deg);
                    [prog_A3, Cell_beta_ca{ioo,ioo_nei}] = sossosvar(prog_A3, beta_deg);
                    [prog_A3, Cell_tau_ca{ioo,ioo_nei}] = sossosvar(prog_A3, tau_deg);
                end
                prog_A3 = sosineq(prog_A3, Cell_Vdec{ioo,1} + kst_var*Lst - L1 );
                prog_A3 = sosineq(prog_A3, s1*(Cell_Vdec{ioo,1}-gamma_in+kst*Lst) - s11*(V0-Eta) );
                prog_A3 = sosineq(prog_A3, (Cell_Vdec{ioo,1} - gamma_in) * s2  - (L2 + Cell_dVdec_dt{ioo,1} ));
                Sto_in = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                        Sto_in = Sto_in + Cell_Vdec{ioo_nei,1}*(Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei})...
                            + Cell_beta_ca{ioo,ioo_nei}*V0 + Cell_tau_ca{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}- gamma_in);
                    else
                        Sto_in = Sto_in + Cell_Vdec{ioo_nei,1}*(Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei})...
                            + Cell_beta_ca{ioo,ioo_nei}*V0 + Cell_tau_ca{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}-gamma_mat_var(ioo_nei,1));
                    end
                end
                prog_A3 = sosineq(prog_A3,  -Cell_dVdec_dt{ioo,1} - L3 + Sto_in );
                prog_A3 = sosineq(prog_A3, (V0 - Eta) * s3_in  - s31_in*(L2 + Vdot0 ));
                ss_AB = 0;
                ss_B = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                    else
                        prog_A3 = sosineq(prog_A3,  (Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei}) - (eps) );
                    end
                    ss_AB = ss_AB + (Cell_alpha_ca{ioo,ioo_nei} - Cell_beta_ca{ioo,ioo_nei});
                    ss_B = ss_B + (Cell_beta_ca{ioo,ioo_nei});
                end
                prog_A3 = sosineq(prog_A3,  -(ss_AB) - (eps) );
                prog_A3 = sosineq(prog_A3,  (ss_B) - (eps) );
                prog_A3 = sossetobj(prog_A3, kst_var);
                prog_A3 = sossolve(prog_A3);
                kst_in = sosgetsol(prog_A3, kst_var);
                kst = kst_in;
                s1_in = sosgetsol(prog_A3,s1);
                s2_in = sosgetsol(prog_A3,s2);
                s11_in = sosgetsol(prog_A3,s11);
                for ioo_nei = 1:1:lf
                    Cell_alpha_in{ioo,ioo_nei} = double(sosgetsol(prog_A3,Cell_alpha_ca{ioo,ioo_nei}));
                    Cell_beta_in{ioo,ioo_nei} = double(sosgetsol(prog_A3,Cell_beta_ca{ioo,ioo_nei}));
                    Cell_tau_in{ioo,ioo_nei} = sosgetsol(prog_A3,Cell_tau_ca{ioo,ioo_nei});
                end            
                prog_A4 = sosprogram(Z);
                p_deg = monomials(Cell_Zdec{ioo,1}, 1);
                [prog_A4, p0_var] = sossosvar(prog_A4, p_deg);
                pdot_z_var = jacobian(p0_var, Z);
                pdot0_var = pdot_z_var * Cell_dZdec_dt{ioo,1};
                prog_A4 = sosineq(prog_A4, p0_var + kst*Lst - L1 );
                prog_A4 = sosineq(prog_A4, s1_in*(p0_var-gamma_mat_var(ioo,1)+kst*Lst) - s11_in*(V0-Eta) );
                prog_A4 = sosineq(prog_A4, (p0_var - gamma_mat_var(ioo,1)) * s2_in  - (L2 + pdot0_var ));
                Sto_in = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                        Sto_in = Sto_in + p0_var*(Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei})...
                            + Cell_beta_in{ioo,ioo_nei}*V0 + Cell_tau_in{ioo,ioo_nei}*(p0_var-gamma_mat_var(ioo_nei,1));
                    else
                        Sto_in = Sto_in + Cell_Vdec{ioo_nei,1}*(Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei})...
                            + Cell_beta_in{ioo,ioo_nei}*V0 + Cell_tau_in{ioo,ioo_nei}*(Cell_Vdec{ioo_nei,1}-gamma_mat_var(ioo_nei,1));
                    end
                end
                prog_A4 = sosineq(prog_A4, -pdot0_var - L3 + Sto_in );
                prog_A4 = sosineq(prog_A4, (V0 - Eta) * s3_in  - s31_in*(L2 + Vdot0 ));                
                ss_AB = 0;
                ss_B = 0;
                for ioo_nei = 1:1:lf
                    if ioo == ioo_nei
                    else
                        prog_A4 = sosineq(prog_A4,  (Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei}) - (eps) );
                    end
                    ss_AB = ss_AB + (Cell_alpha_in{ioo,ioo_nei} - Cell_beta_in{ioo,ioo_nei});
                    ss_B = ss_B + (Cell_beta_in{ioo,ioo_nei});
                end
                prog_A4 = sosineq(prog_A4,  -(ss_AB) - (eps) );
                prog_A4 = sosineq(prog_A4,  (ss_B) - (eps) );
                prog_A4 = sossolve(prog_A4);                
                p0_in = 0.1*sosgetsol(prog_A4, p0_var)+0.9*Cell_Vdec{ioo,1};
                Cell_Vdec{ioo,1} = p0_in/gamma_mat_var(ioo,1);
                Cell_dVdec_dt{ioo,1} = jacobian(Cell_Vdec{ioo,1}, Z) * Cell_dZdec_dt{ioo,1};
                Store_A1_gamma = [Store_A1_gamma; gamma_mat_var(ioo,1)];
                Store_p = [Store_p; p0_in/gamma_mat_var(ioo,1)];
                if i_iter2 > 1
                    p_error = Store_p(end) - Store_p(end-1);
                    p_error_value = double(max(p_error.coefficient+0));
                    Sto_error_value = [Sto_error_value;p_error_value]
                    if double(max(p_error.coefficient+0)) < eps
                        break;
                    end
                end
            end
        end
        % Construct the global Lyapunov function
        A = [];
        B = [];
        for i = 1:1:15
            for j = 1:1:15
                A(i,j) = double(Cell_alpha_in{i,j});
                B(i,j) = double(Cell_beta_in{i,j});
            end
        end
        Out1 = A-B;
        max(real(eig(Out1)));
        sum(Out1,2);
        Out2 = -(A-B);
        max(real(eig(Out2)));
        sum(Out2,2);
        Outt = [sum(Out1,2),sum(Out2,2)];
        out_inv = Out1 - Out1';
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
            Vg = Vg + Cell_Vdec{iop,1}*w_matrix(iop,1);
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
