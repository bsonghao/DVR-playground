# Model Hamiltonian
$$
\begin{aligned}
&\hat{H} = \hat{T} + \hat{V}\\
&\hat{T} = \hat{T}_1\otimes\hat{1} + \hat{1}\otimes\hat{T}_2\\
&\hat{T}_k = \frac{c_k\hat{p}_k^2}{2},\quad k=1,2\\
&\hat{V} = \gamma_1\hat{Q}_1^2\otimes\hat{1} + \epsilon_1\hat{Q}_1^4\otimes\hat{1} + \hat{1}\otimes\kappa_2\hat{Q}_2+\hat{1}\otimes\gamma_2\hat{Q}_2^2+\gamma_{12}\hat{Q}_1\otimes\hat{Q}_2
\end{aligned}
$$
# 1. Contruct the DVR basis
Since the matrix $\bf{Q}$ is in tri-diagonal form, we can construct the proper DVR through diagonalization:

* step 1: Choose $m_k=1$ and $\omega_k=1$ to keep the scaling consistent and numerically stable

* step 2: Choose equilibrium shifts $Q_{0,k}$ based on V:
 
   If we ignore the mode coupling terms
   * For Q1 part, it is a double well potential:
$$
V_1(Q_1) = \gamma_1Q_1^2+\epsilon_1Q_1^4,\quad Q_{1,min}\approx\pm\sqrt{\frac{-\gamma_1}{2\epsilon_1}}
$$
   
   * For Q2 part, it is a displace harmonic potential:
   $$
   V_2(Q_2) = \gamma_2Q_2^2+\kappa_2Q_2, \quad Q_{2,min}\approx -\frac{\kappa_2}{2\gamma_2}
   $$
   
   * Coupling $\gamma_{12}>0$ favors opposite signs for Q1 and Q2 to lower energy, so a reasonable single-well expansion point is
   $$
   Q_{0,1} \approx -3.873, \quad Q_{0,2}\approx 0.625
   $$

* step 3: evalulation matrix element in H.O. basis analytically:
$$
Q_{ij}^{k} = \langle\varphi_i^{(k)}\mid\hat{Q}_k\mid\varphi_j^{(k)}\rangle=Q_0^{(k)}\delta_{ij} + \delta_{i,j-1}\sqrt{\frac{j}{2m_k\omega_k}}+\delta_{i,j+1}\sqrt{\frac{j+1}{2m_k\omega_k}}
$$
Convert to dimensionless Q: $ \quad \tilde{Q}_k = \hat{Q}_k\sqrt{m_k\omega_k}$, we have:
$$
\tilde{Q}_{ij}^{k} = \sqrt{m_k\omega_k}\langle\varphi_i^{(k)}\mid\hat{Q}_k\mid\varphi_j^{(k)}\rangle=\tilde{Q}_0^{(k)}\delta_{ij} + \delta_{i,j-1}\sqrt{\frac{j}{2}}+\delta_{i,j+1}\sqrt{\frac{j+1}{2}}
$$
* step 4: Diagonalize the $\bf{Q}$ matrix, its eigenvalues are the DVR points $\{Q_{\alpha}^{(k)}\}_{\alpha=1}^{N_k}$ and the normalized eigenvectors define the HO$\rightarrow$DVR transform
$$
\bf{\tilde{Q}}^{(k)} = \bf{U}^{(k)}\bf{X}^{(k)}\bf{U}^{(k)\dagger},\quad \bf{X}_{\alpha\beta}^{(k)}=\tilde{Q}_{\alpha}^{(k)}\delta_{\alpha\beta}
$$

$$
\mid\theta^{(k)}\rangle = \bf{U}^{(k)\dagger}\mid\varphi^{(k)}\rangle
$$

* step 5. Obtain DVR grids and the corresponding weights:
 
  * For symmetric tri-diagonal "Jacobi" constructions, the quadrature weights are the squares of the first components of the normalized eigenvectors:
$$
w_{\alpha_k}^{(k)} = \left(\frac{U_{1,\alpha_k}^{(k)}}{\varphi_1(Q_{\alpha_k})}\right)^2
$$
where
$$
\varphi_1(Q_{\alpha_k}) = \pi^{-\frac{1}{4}}\exp\left(-\frac{1}{2}\left(\tilde{Q}_{\alpha_k}^{(k)}-\tilde{Q}_0^{(k)}\right)^2\right)
$$
  * Enforce positive weights by flipping the sign of any eigenvector whose first component is negative to remove the arbitrary phase.
  
**Verify the quadrature rule from the weights and grids obtained from diagonalization DVR**:

The weights and grids are supposed to obey the Gaussian quadratures rule:
\begin{align*}
\langle\varphi_i^{(k)}\mid\tilde{Q}^{(k)}\mid\varphi_j^{(k)}\rangle &= \int w(\tilde{Q}^{(k)})\left(\tilde{Q}^{(k)}\right)^lP^*_{i-1}(\tilde{Q}^{(k)})P_{j-1}(\tilde{Q}^{(k)})d\tilde{Q}^{(k)}\\
&=\sum_{\alpha_k=1}^{N_{k}}w_{\alpha_k}^{(k)}\varphi_i^*(\tilde{Q}_{\alpha_k}^{(k)})\left(\tilde{Q}_{\alpha_k}^{(k)}\right)^l\varphi_j(\tilde{Q}_{\alpha_k}^{(k)})
\end{align*}



Where $\forall \ i-0,1,...,N_k-1$:
$$
\varphi_i(\tilde{Q}_{\alpha_k}^{(k)})=\frac{1}{\sqrt{2^ii!}}\pi^{-\frac{1}{4}}\exp\left(-\frac{1}{2}(\tilde{Q}_{\alpha_k}^{(k)}-\tilde{Q}_0^k)^2\right)H_i\left(\tilde{Q}_{\alpha_k}^{(k)}-\tilde{Q}_0^k\right)
$$
The quadrature rule is exact $\forall \ i,j,l: \ i+j+l\leq 2N_{k}+1$ or for $l=0,1$ the quadrature rule is supposed to be exact.

# 2. Construct the Hamiltonian matrix

* Construct kinetic energy term in DVR
  * In H.O. basis, the H.O. Hamiltonian is diagonal:
  $$
  \hat{H}_{H.O.}^{(k)} = \frac{\hat{p}_k^2}{2m_k} + \frac{1}{2}m_k\omega_k^2\hat{Q}_k^2, \quad \langle\varphi_i^{(k)}\mid\hat{H}_{H.O.}^{(k)}\mid\varphi_j^{(k)}\rangle = \omega_k(j+\frac{1}{2})\delta_{ij}
  $$
  
  * Solve for $\hat{P}_k^2$ in H.O. basis:
  $$
  \hat{P}^{2}_{k} = 2m_k\hat{H}_{H.O.}^{(k)} - m_k^2\omega_k^2\hat{Q}_k^2.
  $$
  Convert int dimension less momentum and position opertors:
  $$
  \tilde{P}_k = \frac{\hat{P}_k}{\sqrt{m_k\omega_k}}, \quad \tilde{Q}_k = \hat{Q}_k\sqrt{m_k\omega_k}
  $$
  We have:
  $$
\bf{\tilde{P}^2}^{(k)}_{H.O.} = \frac{2\bf{H}_{H.O.}^{(k)}}{\omega_k}-\left(\tilde{\bf{Q}}_k^2\right)_{H.O.}
  $$
  * Transform to DVR:
  $$
  \bf{\tilde{P}^2}^{(k)}_{DVR} = \bf{U}^{(k)\dagger}\bf{\tilde{P}^2}^{(k)}_{H.O.}\bf{U}^{(k)}
  $$
  
 * Potential energy in the product DVR is diagonal:
 
On each grid point $(Q_{\alpha_1}^{(1)}, Q_{\alpha_2}^{(2)})$ compute:
   $$
   \bf{V}_{DVR}^{\alpha_1,\alpha_2} = \gamma_1\left(Q_{\alpha_1}^{(1)}\right)^2\otimes\bf{I}_{N_2} + \varepsilon_1\left(Q_{\alpha_1}^{(1)}\right)^4\otimes\bf{I}_{N_2} + \bf{I}_{N_1}\otimes\kappa_2\left(Q^{(2)}_{\alpha_2}\right) + \bf{I}_{N_1}\otimes\gamma_2\left(Q^{(2)}_{\alpha_2}\right)^2 + \gamma_{12}Q_{\alpha_1}^{(1)}\otimes Q_{\alpha_2}^{(2)}.
   $$
   
   * Assemble the full Hamiltonian in theproduct DVR using Kronecker sums / products:
   $$
   \bf{H}_{DVR} = \bf{T}_{DVR}^{(1)}\otimes\bf{I}_{N_2} + \bf{I}_{N_1}\otimes\bf{T}_{DVR}^{(2)} + \bf{V}_{DVR}
   $$
   
* Diagonalize H and report the five lowest eigenvalues.
$$
\bf{H}_{DVR}\bf{C} = \bf{C}\bf{E},\quad \bf{E}_{ij}=\varepsilon_i\delta_{ij}
$$

# 3. Wave functions on the DVR grid

* Using the discrete-$\delta$ property $\theta_{\alpha_k}^{(k)}(Q_{\beta_k}^{(k)})= \frac{\delta_{\alpha_k\beta_k}}{\sqrt{w_{\alpha_k}^{(k)}}}$:
\begin{align*}
&\psi_n(Q_{\alpha_1}^{(1)},Q_{\alpha_2}^{(2)})=\sum_{\beta_1,\beta_2}C^n_{\beta_1, \beta_2}\theta_{\beta_1}^{(1)}(Q_{\alpha_1}^{(1)})\theta_{\beta_2}^{(2)}(Q_{\alpha_2}^{(2)})\\
&=\frac{\sum_{\beta_1,\beta_2}C^n_{\beta_1,\beta_2}\delta_{\beta_1,\alpha_1}\delta_{\beta_2,\alpha_2}}{\sqrt{w_{\beta_1}^{(1)}}\sqrt{w_{\beta_2}^{(2)}}}\\
&=\frac{C^n_{\alpha_1,\alpha_2}}{\sqrt{w_{\alpha_1}^{(1)}}\sqrt{w_{\alpha_2}^{(2)}}},
\end{align*}
where
$$
C_{\alpha_1,\alpha_2}^n = \langle\theta_{\alpha_1}^{(1)}\theta_{\alpha_2}^{(2)}\mid\psi_n\rangle 
$$
is the reshaped eigenvector of $\bf{H}_{DVR}$.
