clear; close all;

addpath(genpath("./utils"));

input_info.if_save_data = true;
input_info.if_save_iteration = false;
input_info.sig0 = 1e2;
initial_type = "rand";

prefix1 = "BQP-r1/20-1";
prefix = "admmdyn-data" + "/" + prefix1 + "/";

%% load data
SDP_data = load(prefix + "SDP_data");
At = SDP_data.At;
b = SDP_data.b;
C = SDP_data.C;
input_info.dirname = prefix;
input_info.if_sigfix = true;

%% generate input_info.SDP_data
A = At';
[R, ~, P] = chol(A * At);
SDP_data.A = A; SDP_data.At = At; SDP_data.b = b;
SDP_data.C = C; SDP_data.R = R; SDP_data.P = P;
input_info.SDP_data = SDP_data;
vec_len = size(At, 1);
mat_size = (-1 + sqrt(1 + 8 * vec_len)) / 2;
SDP_data.vec_len = vec_len;
SDP_data.mat_size = mat_size;

%% generate initial X, y, S
if initial_type == "zero"
    input_info.X0 = zeros(size(C));
    input_info.y0 = zeros(size(b));
    input_info.S0 = zeros(size(C));
elseif initial_type == "rand"
    input_info.X0 = randn(size(C));
    input_info.y0 = randn(size(b));
    input_info.S0 = randn(size(C));
elseif initial_type == "warmstart"
    load(prefix + "result");
    input_info.X0 = full(svec_single(data.X_mat));
    input_info.y0 = data.y;
    input_info.S0 = full(svec_single(data.S_mat));
end

%% generate input_info.ADMM_info
ADMM_info.maxiter = 100e4;
ADMM_info.tol = 1e-10;
ADMM_info.scaleA = 0;
ADMM_info.scaleData = 0;
ADMM_info.sigfix = input_info.if_sigfix;
ADMM_info.rank_threshold = 1e-10;
ADMM_info.ADMM_warmstart_iter = ADMM_info.maxiter;
ADMM_info.save_loop = 1e4;
input_info.ADMM_info = ADMM_info;

%% call ADMM
output_info = ADMM_single(input_info);
X_mat = smat_single(output_info.X);
S_mat = smat_single(output_info.S);
Xb_mat = X_mat - output_info.data.sig * S_mat;
pinf_list = output_info.data.pinf_list;
dinf_list = output_info.data.dinf_list;
relgap_list = output_info.data.relgap_list;
Xb_rank_list = output_info.data.Xb_rank_list;
Xb_diff_norm_next_list = output_info.data.Xb_diff_norm_next_list;

%% plot
options.rows = 350;
options.cols = 350;
options.linewidth = 1.5;
options.ms = 12;
options.xls = 13;
options.yls = 13;
options.fontsize = 13;
linecolors = linspecer(2, 'qualitative');
LineColors = flipud(linecolors);
options.LineColors = LineColors;
draw_rmax_dZ_with_mineig(output_info.data, options);

%% helper functions
function output_info = ADMM_single(input_info)
    % X, S are given as svec
    X = input_info.X0;
    y = input_info.y0;
    S = input_info.S0;
    sig = input_info.sig0;
    
    if_save_data = input_info.if_save_data;
    if_save_iteration = input_info.if_save_iteration;
    dirname = input_info.dirname;
    if if_save_data
        dirname = input_info.dirname;
        if ~exist(dirname)
            mkdir(dirname);
        end
    end
    
    ADMM_info = input_info.ADMM_info;

    SDP_data = input_info.SDP_data;
    A = SDP_data.A; At = SDP_data.At; b = SDP_data.b;
    C = SDP_data.C; P = SDP_data.P; R = SDP_data.R; 
    borg = b;
    Corg = C;
    normborg = 1 + norm(borg);
    normCorg = 1 + norm(Corg);
    
    % scale At
    if ADMM_info.scaleA
        normA = sum(At .* At, 1)';
        normA = max(1.0, sqrt(normA));
        DA = spdiags(1./normA, 0, size(At, 2), size(At, 2));
        At = At * DA;
        A = At';
        b = b ./ normA;
        y = normA .* y;
        [R, ~, P] = chol(A * At);
    else
        normA = 1.0;
    end

    % scale b, C
    normb = 1 + norm(b);
    normC = 1 + norm(C);
    if ADMM_info.scaleData
        bscale = normb;
        Cscale = normC;
    else
        bscale = 1;
        Cscale = 1;
    end
    objscale = bscale * Cscale;
    b = b / bscale;
    C = C / Cscale;
    X = X / bscale;
    S = S / Cscale;
    y = y / Cscale;

    if if_save_data
        saved_SDP_data.A = A; saved_SDP_data.At = At; saved_SDP_data.b = b;
        saved_SDP_data.C = C; saved_SDP_data.R = R; saved_SDP_data.P = P;
        saved_SDP_data.normA = normA; 
        saved_SDP_data.bscale = bscale; saved_SDP_data.Cscale = Cscale;
        save(dirname + "normalized_SDP_data.mat", "saved_SDP_data");
    end
    
    % main loop
    prim_win = 0;
    dual_win = 0;
    rescale = 1;
    maxfeas = Inf;
    pinf_list = [];
    dinf_list = [];
    relgap_list = [];
    pobj_list = [];
    dobj_list = [];

    % two steps ahead to get Xb_prev and Xb
    rhsy = 1/sig * b - A * (1/sig * X + S - C);
    y = chol_solve(rhsy, R, P);
    Xb_prev = X + sig * (At * y - C); 
    Xb_proj_prev = proj_single(Xb_prev);
    S = 1/sig * (Xb_proj_prev - Xb_prev);
    X = X + sig * (S + At * y - C);

    rhsy = 1/sig * b - A * (1/sig * X + S - C);
    y = chol_solve(rhsy, R, P);
    Xb = X + sig * (At * y - C); 
    Xb_proj = proj_single(Xb);
    S = 1/sig * (Xb_proj - Xb);
    X = X + sig * (S + At * y - C);

    % three buffers
    Xb_diff_norm_next_buffer = [];
    Xb_diff_ang_triple_buffer = [];
    Xb_rank_buffer = [];
    Xb_diff_norm_next_list = [];
    Xb_diff_ang_triple_list = [];
    Xb_rank_list = [];
    if_restart = false;

    for iter = 1: ADMM_info.maxiter
        % rescale the whole problem
        normy = norm(y);
        normAty = norm(At * y);
        normX = norm(X);
        normS = norm(S);
        normyS = max([normy, normAty, normS]);

        rhsy = 1/sig * b - A * (1/sig * X + S - C);
        y = chol_solve(rhsy, R, P);
        Xb_next = X + sig * (At * y - C);
        [Xb_proj_next, d] = proj_single(Xb_next);
        S = 1/sig * (Xb_proj_next - Xb_next);
        X = X + sig * (S + At * y - C);
        
        % update three buffers
        Xb_diff = Xb - Xb_prev; Xb_diff_next = Xb_next - Xb;
        Xb_diff_norm = norm(Xb_diff_next);
        cos_ang = (Xb_diff' * Xb_diff_next) / (norm(Xb_diff) * norm(Xb_diff_next));
        cos_ang = min(cos_ang, 1.0 - 1e-16);
        cos_ang = max(cos_ang, -1.0 + 1e-16);
        ang = acos(cos_ang);
        d_pos = d(d > 0);
        r = nnz(d_pos > ADMM_info.rank_threshold);
        Xb_prev = Xb;
        Xb = Xb_next;

        Xb_diff_norm_next_list =  [Xb_diff_norm_next_list, Xb_diff_norm];
        Xb_diff_ang_triple_list = [Xb_diff_ang_triple_list, ang];
        Xb_rank_list = [Xb_rank_list, r];

        pinf = norm(normA .* (A * X - b)) * bscale / normborg;
        dinf = norm(At * y + S - C) * Cscale / normCorg;
        pobj = full(C' * X * objscale);
        dobj = b' * y * objscale;
        relgap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));

        fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                iter, Xb_diff_norm, ang, r, pinf, dinf, relgap, pobj, dobj, sig);
        maxfeas = max([pinf, dinf]);
        pinf_list = [pinf_list, pinf];
        dinf_list = [dinf_list, dinf];
        relgap_list = [relgap_list, relgap];
        pobj_list = [pobj_list, pobj];
        dobj_list = [dobj_list, dobj];
        if max([pinf, dinf, relgap]) < ADMM_info.tol
            fprintf("Convergent! \n");
            break;
        end

        % save each iteration's Xb 
        if if_save_iteration
            data.X_mat = smat_single(bscale * X);
            data.y = Cscale * (y ./ normA);
            data.S_mat = smat_single(Cscale * S);
            data.sig = sig;
            data.Xb_mat = data.X_mat - sig * data.S_mat;
            Xb_mat_iter = data.Xb_mat;
            save(dirname + sprintf("%d.mat", iter), "Xb_mat_iter");
        end

        % tmp save
        if rem(iter, ADMM_info.save_loop) == (ADMM_info.save_loop - 1)
            filename = dirname + "result";
            data.pinf_list = pinf_list;
            data.dinf_list = dinf_list;
            data.relgap_list = relgap_list;
            data.pobj_list = pobj_list;
            data.dobj_list = dobj_list;
            data.Xb_diff_norm_next_list = Xb_diff_norm_next_list;
            data.Xb_diff_ang_triple_list = Xb_diff_ang_triple_list;
            data.Xb_rank_list = Xb_rank_list;

            data.X_mat = smat_single(bscale * X);
            data.y = Cscale * (y ./ normA);
            data.S_mat = smat_single(Cscale * S);
            data.sig = sig;
            data.Xb_mat = data.X_mat - sig * data.S_mat;

            data.min_sigular_val_Xb = min(abs(eig(data.Xb_mat)));
            data.ADMM_info = ADMM_info;
            if if_save_data
                save(filename, "data");
            end
        end
    end

    output_info.X = bscale * X;
    output_info.y = Cscale * (y ./ normA);
    output_info.S = Cscale * S;
    output_info.pobj = pobj;
    output_info.dobj = dobj;
    X = output_info.X;
    y = output_info.y;
    S = output_info.S;

    filename = dirname + "result";
    data.pinf_list = pinf_list;
    data.dinf_list = dinf_list;
    data.relgap_list = relgap_list;
    data.pobj_list = pobj_list;
    data.dobj_list = dobj_list;
    data.Xb_diff_norm_next_list = Xb_diff_norm_next_list;
    data.Xb_rank_list = Xb_rank_list;

    data.X_mat = smat_single(X);
    data.y = y;
    data.S_mat = smat_single(S);
    data.sig = sig;
    data.Xb_mat = smat_single(X - sig * S);

    data.min_sigular_val_Xb = min(abs(eig(data.Xb_mat)));
    data.ADMM_info = ADMM_info;
    if if_save_data
        save(filename, "data");
    end

    output_info.data = data;
end

function y = chol_solve(rhsy, R, P)
    rhsy = P' * rhsy;
    tmp = R' \ rhsy;
    tmp = R \ tmp;
    y = P * tmp;
end

function draw_rmax_dZ_with_mineig(data, options)
    rows = options.rows;
    cols = options.cols;
    linewidth = options.linewidth;
    ms = options.ms;
    xls = options.xls;
    yls = options.yls;
    fontsize = options.fontsize;
    LineColors = options.LineColors;

    figure('Position', [0, 0, cols, rows]);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    relgap_list = data.relgap_list;
    pinf_list = data.pinf_list;
    dinf_list = data.dinf_list;
    dZ_list = data.Xb_diff_norm_next_list;
    rmax_list = max([relgap_list; pinf_list; dinf_list]);


    semilogy(rmax_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    semilogy(dZ_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on;
    grid on;
    legend({'$r_{\max}^{(k)}$', '$\parallel Z^{(k+1)} - Z^{(k)} \parallel_F$'}, 'Interpreter', 'latex', 'FontSize', fontsize);


    mineig = data.min_sigular_val_Xb;
    exponent = floor(log10(abs(mineig)));
    mantissa = mineig / 10^exponent;
    tmp = sprintf('$\\lambda_{\\min}(|Z_{\\star}|) = %.1f \\times 10^{%d}$', mantissa, exponent);
    title(tmp, 'FontSize', xls, 'Interpreter', 'latex');
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');

    axis([-inf inf -inf inf]);
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = xls; % Set the x-axis tick font size
    ax.YAxis.FontSize = yls; % Set the x-axis tick font size
    ax.GridLineWidth = 0.6; % Set the grid line width 
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.4;
end



















