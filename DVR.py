import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import itertools as it
from scipy.special import hermite, factorial

def HO_basis(n, x, qe):
    """compute H.O. basis function"""
    H_n = hermite(n-1, monic=False)
    delta = x-qe
    phi = H_n(delta) * np.exp(-0.5*delta**2)
    phi /= (np.sqrt(2.0**(n-1)*factorial(n-1))*np.pi**0.25)
    return phi

class DVR_problem(object):
    """
    solve the DVR problem set
    """

    def __init__(self, c1, c2, g1, eps1, k2, g2, g12):
        """
        c1,c2,gamma1, eps1, k2, g2, g12 are hamiltonian parameters
        """

        self.c1=c1
        self.c2=c2
        self.g1 = g1
        self.eps1 = eps1
        self.k2 = k2
        self.g2 = g2
        self.g12 = g12


    def construct_DVR(self, omega, m, N):
        """
        construct DVR from H.O. basis

        input
        ---------
        omega
           a list that stores the frequencies of the two modes of the H.O. basis
        m: list
           a list that stores the mass of the two modes for the H.O. basis.
        N: list
           a list that stores the truncation of the H.O. basis for the two modes

        output
        ------
        DVR: a dictionary that store DVR with key of grid points and values of weights
        """

        def make_Q_matrix(n, q_e):
            """
            construct Q matrix
            """
            Q = np.zeros((n, n), dtype=float)
            Q += np.eye(n) * q_e
            # Off diagonal elements
            off = np.sqrt(np.arange(1, n) / 2)
            Q[np.arange(n-1), np.arange(1, n)] = off # upper
            Q[np.arange(1, n), np.arange(n-1)] = off # lower

            #  Q_ = Q.copy()
            #  for i,j in it.product(range(n), repeat=2):
            #      if i == j+1:
            #          Q_[i, j] = np.sqrt((j+1)/(2*mw))
            #      elif i == j-1:
            #          Q_[i, j] = np.sqrt(j/(2*mw))
            #      else:
            #          pass

            #  assert np.allclose(Q, Q_)
            #  assert np.allclose(Q_, Q_.transpose())
            assert np.allclose(Q, Q.transpose())

            return Q



        # determine the equlibrium position
        Q_0 = np.zeros(2)
        Q_0[0] = -np.sqrt(-self.g1 / (2*self.eps1))
        Q_0[1] = -self.k2 / (2*self.g2)

        print(f"Equlibrium position for the two H.O. basis:\n{Q_0}")

        self.Qpts = {}
        self.Wts = {}
        self.U = {}
        self.Q_HO = {}
        # loop over to two modes
        for mode in range(2):
            # construct the Q matrix
            mw = m[mode]*omega[mode]
            Q_HO = make_Q_matrix(N[mode], Q_0[mode])

            # store Q matrix
            self.Q_HO[mode] = Q_HO.copy()

            # diagonalize Q matrix
            evals, evecs = la.eigh(Q_HO)

            # the eigenvalues are the DVR grid points
            self.Qpts[mode] = evals.copy()

            # quadrature weights = squares of first components of normalized eigenvector
            vecs = evecs.copy()
            signs = np.sign(vecs[0, :])
            signs[signs == 0] = 1.0
            vecs *= signs  # flip columns where needed
            # print(vecs[0,:])
            psi_1 = HO_basis(1, self.Qpts[mode], Q_0[mode])
            if True:
                # check if the psi_1 is calculated properly
                assert np.allclose(psi_1, np.exp(-0.5*(self.Qpts[mode]-Q_0[mode])**2)*(1./np.pi)**0.25)
            # print(psi_1)
            machine_precision = 1e-9
            self.Wts[mode] = ((vecs[0,:])/(psi_1))**2

            # store transform: HO->DVR
            self.U[mode] = vecs.copy()

            print(f"For mode {mode+1}:")
            print(f"DVR grids: \n{self.Qpts[mode]}")
            print(f"DVR weights: \n{self.Wts[mode]}")


            # check the quadrature rule
            if True:
                # evaluate matrix element from quadrature rule
                delta = np.zeros(shape=(N[mode], N[mode]), dtype=float)
                Q_quad = np.zeros(shape=(N[mode], N[mode]), dtype=float)
                for i,j in it.product(range(N[mode]), repeat=2):
                    phi_i = HO_basis(i+1, self.Qpts[mode], Q_0[mode])
                    phi_j = HO_basis(j+1, self.Qpts[mode], Q_0[mode])
                    delta[i, j] = sum(self.Wts[mode]*phi_i*phi_j)
                    Q_quad[i, j] = sum(self.Wts[mode]*phi_i*self.Qpts[mode]*phi_j)

                # check if the matrix elements obey exact quadrature rule for i+j+l<2n+1
                assert np.allclose(delta, np.eye(N[mode])) # for l=0 <phi_i|phi_j>_{quad} = delta_ij
                assert np.allclose(Q_quad, Q_HO) # for l=1 <phi_i|x|phi_j>_{quad} = X_{HO}

        return

    def construct_H(self, m, omega, N, lowest=5):
        """
        Construct Hamitonian matrix in DVR
        """
        N1, N2 = N[0], N[1]
        c_list = [self.c1, self.c2]
        g1, eps1, k2, g2, g12 = self.g1, self.eps1, self.k2, self.g2, self.g12

        K_DVR = {}

        # firstly construct KEO in DVR
        for mode in range(2):

            # unpack parameters
            Q_HO = self.Q_HO[mode].copy()
            U = self.U[mode]
            c = c_list[mode]

            # Obtain P2_HO form analytical expression of HO Hamiltonian in H.O. basis
            H_HO = np.diag((np.arange(N[mode])+0.5))
            P2_HO = 2*H_HO.copy()
            P2_HO -= (Q_HO @ Q_HO)
            # transfrom from HO to DVR basis
            P2_DVR = U.T @ P2_HO @ U
            K_DVR[mode] = 0.5 * c * P2_DVR.copy()

        # full KEO is is obtained via Kronecker sums
        I1 = np.eye(N[0])
        I2 = np.eye(N[1])
        K = np.kron(K_DVR[0], I2) + np.kron(I1, K_DVR[1])
        print(f"Full KEO: {K.shape}")

        # secondly construct potential operator

        Q1, Q2 = np.meshgrid(self.Qpts[0], self.Qpts[1], indexing="ij")
        V = g1 * Q1**2 + eps1*Q1**4 + k2*Q2 + g2*Q2**2 + g12 * Q1 * Q2
        print(f"Full potential operator: {V.shape}")
        self.V = V.copy()

        # check if the potential term is constructed properly
        if True:
            Q1, Q2 = np.diag(self.Qpts[0]), np.diag(self.Qpts[1])
            V_ = np.kron(g1*Q1@Q1 + eps1*Q1@Q1@Q1@Q1, I2) + np.kron(I1, k2*Q2 + g2*Q2@Q2) + g12*np.kron(Q1, Q2)
            assert np.allclose(np.diag(V.ravel()), V_)

        # The Full Hamiltonian= K + V
        self.H_DVR = K + np.diag(V.ravel())
        print(f"Full Hamiltonian: {self.H_DVR.shape}")
        # check if the H is symmetric
        assert np.allclose(self.H_DVR, self.H_DVR.transpose())

        E, C = la.eigh(self.H_DVR)

        # sort eigenvalues and eigenvectors
        indices = np.argsort(E)
        E = E[indices]
        C = C[:,indices]

        # store eigenvectors
        self.C = C

        # report first few lowest eigenvalues
        E5 = E[:lowest]

        for i, e in enumerate(E5, 1):
            print(f"{i}: {e:.8f}")

        return

    def plot_DVR_wfn(self, m, omega, N, lowest=1):
        """
        plot eigen function in DVR representation
        """

        # obtain wfn for first few lowest states from DVR basis
        machine_precision = 1e-9
        sqrtw1 = 1./(np.sqrt(self.Wts[0])+machine_precision)  # (N1,1)
        sqrtw2 = 1./(np.sqrt(self.Wts[1])+machine_precision)  # (1,N2)
        psi = []
        for n in range(lowest):
            coeff = self.C[:,n].reshape(N[0], N[1])
            # print(coeff)
            # print(sqrtw1)
            # print(sqrtw2)
            wfn = np.einsum("ij,i,j->ij", coeff , sqrtw1 , sqrtw2)
            # print(wfn)
            psi.append(wfn)

        # Plot the first five wavefunctions
        fig, axs = plt.subplots(1, lowest, figsize=(40*lowest, 40), constrained_layout=True)
        levels = 80
        for n in range(lowest):
            ax = axs[n] if lowest > 1 else axs
            cn = ax.contourf(self.Qpts[0], self.Qpts[1], psi[n].T, levels=levels, cmap="coolwarm")
            ax.set_title(f"$\psi_{n+1}(Q_1,Q_2)$", fontsize=200)
            ax.set_xlabel("Q1", fontsize=200)
            ax.set_ylabel("Q2", fontsize=200)
            # ax.set_xlim(-6,-3)
            # ax.set_ylim(9, 11.24)
            ax.tick_params(labelsize=80)
            cbar = fig.colorbar(cn, ax=ax)
            cbar.ax.tick_params(labelsize=160)
        plt.savefig("Plot_DVR_wfn.png")
        plt.show()

        # Plot the potential in DVR representation

        fig, axs = plt.subplots(1, 1, figsize=(40, 40), constrained_layout=True)
        levels = 40

        cn = axs.contourf(self.Qpts[0], self.Qpts[1], self.V.T, levels=levels, cmap="coolwarm")
        axs.set_title(f"$V(Q_1,Q_2)$", fontsize=200)
        axs.set_xlabel("Q1", fontsize=200)
        axs.set_ylabel("Q2", fontsize=200)
        axs.tick_params(labelsize=80)
        #axs.set_xlim(-6,-3)
        #axs.set_ylim(9, 11)
        cbar = fig.colorbar(cn, ax=axs)
        cbar.ax.tick_params(labelsize=80)
        plt.savefig("Plot_DVR_potential.png")
        plt.show()

        return
