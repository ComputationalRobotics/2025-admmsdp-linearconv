clear; close all;

addpath(genpath("./utils"));

options.rows = 250;
options.cols = 250;
options.linewidth = 1.5;
options.ms = 12;
options.xls = 13;
options.yls = 13;
options.fontsize = 13;
linecolors = linspecer(2, 'qualitative');
LineColors = flipud(linecolors);
options.LineColors = LineColors;

prefix1_list = ["BQP-r1/20-1"];

% prefix1_list = [
%     "BQP-r1-zero/30-3";
%     "quasar/500-1";
%     "XM/15";
%     "toy";
% ];

% BQP_1_list = [10, 20, 30, 40];
% BQP_2_list = [1, 2, 3];
% prefix1_list = [];
% for i = 1: length(BQP_1_list)
%     for j = 1: length(BQP_2_list)
%         BQP_1 = BQP_1_list(i);
%         BQP_2 = BQP_2_list(j);
%         prefix1_list = [prefix1_list; sprintf("BQP-r2/%d-%d", BQP_1, BQP_2)];
%     end
% end

% prefix1_list = [
%     "1dc-1024";
%     "cancer-100";
%     "cnhil8";
%     "cnhil10";
%     "G40-mb";
%     "hand";
%     "neosfbr25";
%     "rose13";
%     "swissroll";
% ];

% prefix1_list = [
%     "maxcut/G1";
%     "maxcut/G5";
%     "maxcut/G9";
%     % "maxcut/G11";
%     "maxcut/G14";
%     "maxcut/G18";
% ];

% prefix1_list = [
%     "quasar/100-1";
%     "quasar/200-1";
%     "quasar/500-1";
% ];

% prefix1_list = [
%     "XM/48";
%     "XM/93";
%     "XM/149";
% ];

% prefix1_list = [
%     "hamming-9-8";
%     "hamming-10-2";
%     "hamming-11-2";
% ];

% prefix1_list = [
%     "theta-12";
%     "theta-102";
%     "theta-123";
% ];

for i = 1: length(prefix1_list)
    prefix1 = prefix1_list(i);
    load("./admmdyn-data/" + prefix1 + "/result.mat");

    img_dir = "./figs/" + prefix1;

    if ~exist(img_dir, "dir")
        mkdir(img_dir);
    end

    % options.rows = 100;
    % options.cols = 250;
    % xline_iter = draw_rank(data, img_dir + "/X_rank.png", options, true);
    % options.rows = 150;
    % options.cols = 310;
    % draw_dZ(data, img_dir + "/dZ.png", options, xline_iter);

    img_name = img_dir + "/rmax_dZ.png";
    draw_rmax_dZ_with_mineig(data, img_name, options);
    % draw_rmax_dZ(data, img_name, options);
    
    % if i == 1
    %     img_name = img_dir + "/legend_rank_dZ.png";
    %     draw_legend(img_name, options);
    % end
end

function xline_iter = draw_rank(data, img_name, options, if_identify)
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

    rank_list = data.Xb_rank_list;
    rank_list(rank_list == 0) = 1;
    semilogy(rank_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    grid on;
    % xlabel('Iteration', 'FontSize', xls, 'Interpreter', 'latex');
    % ylabel('rank $(X^{(k)})$', 'FontSize', xls, 'Interpreter', 'latex');

    if if_identify
        final_rank = rank_list(end);
        i = length(rank_list);
        while rank_list(i) == final_rank 
            i = i - 1;
        end
        xline_iter = i;
        xline(xline_iter, ':', 'Color', [0, 0.5, 0], 'LineWidth', 2, 'HandleVisibility', 'off');
    else 
        xline_iter = -1;
    end

    mineig = data.min_sigular_val_Xb;
    exponent = floor(log10(abs(mineig)));
    mantissa = mineig / 10^exponent;
    tmp = sprintf('$\\lambda_{\\min}(|Z_{\\star}|) = %.1f \\times 10^{%d}$', mantissa, exponent);
    title(tmp, 'FontSize', xls, 'Interpreter', 'latex');

    axis([-inf inf -inf inf]);
    ax = gca; % Get the current axes
    % ax.XAxis.FontSize = xls; % Set the x-axis tick font size
    ax.XTickLabel = [];
    ax.YAxis.FontSize = yls; % Set the x-axis tick font size
    ax.GridLineWidth = 0.6; % Set the grid line width 
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.4;
    print(img_name, '-dpng', ['-r', num2str(300)]);
    close all;
end

function draw_dZ(data, img_name, options, xline_iter)
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

    dZ_list = data.Xb_diff_norm_next_list;
    semilogy(dZ_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on;
    grid on;
    xlabel('Iteration', 'FontSize', xls, 'Interpreter', 'latex');
    % ylabel('$\parallel Z^{(k+1)} - Z^{(k)} \parallel_F$', 'FontSize', xls, 'Interpreter', 'latex');

    if nargin > 3
        xline(xline_iter, ':', 'Color', [0, 0.5, 0], 'LineWidth', 2, 'HandleVisibility', 'off');
    end

    axis([-inf inf -inf inf]);
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = xls; % Set the x-axis tick font size
    ax.YAxis.FontSize = yls; % Set the x-axis tick font size
    ax.GridLineWidth = 0.6; % Set the grid line width 
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.4;
    print(img_name, '-dpng', ['-r', num2str(300)]);
    close all;
end

function draw_rmax_dZ_with_mineig(data, img_name, options)
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
    % legend({'$r_{\max}^{(k)}$', '$\parallel Z^{(k+1)} - Z^{(k)} \parallel_F$'}, 'Interpreter', 'latex', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', fontsize);


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
    print(img_name, '-dpng', ['-r', num2str(300)]);
    close all;
end


function draw_rmax_dZ(data, img_name, options)
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
    legend({'$r_{\max}^{(k)}$ ', '$\parallel Z^{(k+1)} - Z^{(k)} \parallel_F$'}, 'Interpreter', 'latex', 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', fontsize);
    xlabel('Iteration', 'FontSize', xls, 'Interpreter', 'latex');
    axis([-inf inf -inf inf]);
    ax = gca; % Get the current axes
    ax.XAxis.FontSize = xls; % Set the x-axis tick font size
    ax.YAxis.FontSize = yls; % Set the x-axis tick font size
    ax.GridLineWidth = 0.6; % Set the grid line width 
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.4;
    print(img_name, '-dpng', ['-r', num2str(300)]);
    close all;
end

function draw_legend(img_name, options)
    rows = 50;
    cols = 300;
    linewidth = 3;
    ms = options.ms;
    xls = options.xls;
    yls = options.yls;
    fontsize = options.fontsize;
    LineColors = options.LineColors;

    figure('Position', [0, 0, cols, rows]);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    hold on;
    plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    hold off;

    legend({'rank $(X^{(k)})$', '$\parallel Z^{(k+1)} - Z^{(k)} \parallel_F$'}, 'Interpreter', 'latex', 'NumColumns', 2, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off');
    
    set(gca, 'LooseInset', [0,0,0,0]) % Remove extra margins
    set(gca, 'Box', 'off') % Remove box around the axes
    set(gcf, 'Color', 'w') % Set background color to white
    axis off % Hide axes

    print(img_name, '-dpng', ['-r', num2str(300)]);
    close all;
end



