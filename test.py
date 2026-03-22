from DVR import HO_basis, DVR_problem
import numpy as np

def main():
    # ----------------------------
    # Parameters from the problem
    # ----------------------------
    c1, c2   = 0.120, 0.180    # c1, c2
    g1, eps1 = -0.150,0.005   # gamma1, epsilon1
    k2, g2   = -0.010, 0.008   # kappa2, gamma2
    g12      = 0.010           # gamma12

    # Choose basis set truncation HO masses and frequencies (mass-/freq-scaled)
    N = [53, 53]
    m = [2., 8.]
    omega = [0.5, 0.25]

    problem = DVR_problem(c1, c2, g1, eps1, k2, g2, g12)
    problem.construct_DVR(omega, m, N)
    # store DVR grids and weight
    np.savetxt("Q1_points.txt", problem.Qpts[0])
    np.savetxt("W1_weights.txt", problem.Wts[0])
    np.savetxt("Q2_points.txt", problem.Qpts[1])
    np.savetxt("W2_weights.txt", problem.Wts[1])
    problem.construct_H(omega, m, N,lowest=5)
    problem.plot_DVR_wfn(omega, m, N,lowest=5)

    return

if __name__ == '__main__':
    main()
