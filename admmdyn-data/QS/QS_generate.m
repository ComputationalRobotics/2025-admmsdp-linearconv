clear; close all;

addpath("../../pathinfo/");
ksc;
% jordan;
pathinfo("lib") = "../../lib";
pathinfo("sos-sdp-conversion") = "../../sos-sdp-conversion";
pathinfo("utils") = "../../utils";
pathinfo("admmsdp") = "../../admmSDP";
pathinfo("kscutils") = "../../ksc-utils";
pathinfo("modules") = "../../modules";

keys = pathinfo.keys;
for i = 1: length(keys)
    key = keys(i);
    addpath(genpath(pathinfo(key)));
end

n = 40;
id = 1;
kappa = 2;
prefix = sprintf("%d-%d", n, id);
[At, b, C, SDP] = QS(n, kappa);
if ~exist(prefix, "dir")
    mkdir(prefix);
end
save(prefix + "/SDP_data.mat", "At", "b", "C");

%% test Mosek
% prob = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
% [~, res] = mosekopt('minimize info', prob);
% [Xopt, yopt, Sopt, obj] = recover_mosek_sol_blk(res, SDP.sdpt3.blk);
% X_mat_mosek = Xopt{1};


% quartic over sphere
function [At, b, C, SDP] = QS(n, kappa)
    d       = n; % BQP with d variables (when d=20 it is already pretty slow)
    x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
    problem.vars = {[x]};
    % optimize a quartic over a sphere
    tmp = mymonomial(x, 4);
    c = randn(size(tmp, 1), 1);
    problem.objective = c' * tmp;
    problem.inequality = {[msspoly()]};
    problem.equality = {[1 - x' * x]};
    problem.relaxation_order = kappa;
    problem.rip_predecessor = [0];
    problem.regularization.expression = msspoly(1);
    problem.regularization.value = 1;
    [SDP, info] = sparse_sdp_relax(problem);
    
    pool = generate_pool_fast(SDP, 1e-10, 5e-4);
    SDP = remove_redundant(SDP, pool);
    fprintf("\n At's smallest singular value: %3.2e \n", svds(SDP.sedumi.At, 1, 'smallest'));
    At = SDP.sdpt3.At{1}; 
    b = SDP.sdpt3.b; 
    C = svec_single(SDP.sdpt3.C{1}); 
end