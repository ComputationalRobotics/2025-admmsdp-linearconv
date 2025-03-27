## Codes and Data for the Paper [*"Local Linear Convergence of the Alternating Direction Method of Multipliers for Semidefinite Programming under Strict Complementarity"*](https://arxiv.org/abs/2503.20142)

### Installation

There is no need to install anything except MATLAB. 

### Usage 

Please directly run `admmdyn.m` in MATLAB to test `BQP/20-1` data with $\sigma = 100$ and random (standard Gaussian) initial guess. In the end, you will get an image like this:

<p align="center">
  <img src="./admmdyn-data/BQP-r1/20-1/demo.png" width="50%">
</p>

To change to different SDP instances, please change the `prefix1` in line 10, `admmdyn.m` to other tags, such as `theta-12`. The program will automatically load `SDP_data.mat` for that tag. The tag names correspond to labels in paper's Table 1. 

Other important hyper-parameters in `admmdyn.m`:

- `input_info.if_save_data`: whether output information from ADMM optimization procedure will be saved or not. Default: `true`. 
- `input_info.if_save_iteration`: whether save each iteration's $Z^{(k)}$ or not. Default: `false`. Warning: if set to `true`, it will consume large amount of disk storage for even for small-scale problems.
- `input_info.sig0`: fixed $\sigma$'s value in ADMM. Default: $100$. 
- `initial_type`: use all-zero or random initial guess for $(X^{(0)}, y^{(0)}, S^{(0)})$. Default: `"rand"`. 
- `ADMM_info.maxiter`: ADMM iteration limit. Default: $10^6$.
- `ADMM_info.tol`: ADMM max KKT residual tolerance. Default: $10^{-10}$.
- `ADMM_info.scaleA`: whether scale `A` in SDP data. Default: `false`.
- `ADMM_info.scaleData`: whether scale `b` and `C` in SDP data. Default: `false`.

If you have chosen to save the output data, there will be a new `result.mat` file in the `./admmdyn-data` + tag name folder. Its components are: 

- `pinf_list`: a list of $r_p^{(k)}$ during ADMM procedure.
- `dinf_list`: a list of $r_d^{(k)}$ during ADMM procedure. 
- `relgap_list`: a list of $r_g^{(k)}$ during ADMM procedure.
- `pobj_list`: a list of $\langle C, X^{(k)} \rangle$ during ADMM procedure.
- `dobj_list`: a list of $b^T y^{(k)}$ during ADMM procedure. 
- `Xb_diff_norm_next_list`: a list of $|| Z^{(k+1)} - Z^{(k)} ||_F$ during ADMM procedure.
- `Xb_diff_ang_triple_list`: a list of angle between $Z^{(k+1)} - Z^{(k)}$ and $Z^{(k)} - Z^{(k-1)}$ during ADMM procedure. 
- `Xb_rank_list`: a list of $\text{rank} (X^{k})$ during ADMM procedure. 
- `min_sigular_val_Xb`: $\lambda_\min (|Z_\star|)$. 
- `X_mat`, `y`, `S_mat`, `Xb_mat`: converged optimal solutions $X_\star, y_\star, S_\star, Z_\star$.

We also open source the plotting files in `admmdyn_plot.m` for reader to reproduce the plots in the paper. 
