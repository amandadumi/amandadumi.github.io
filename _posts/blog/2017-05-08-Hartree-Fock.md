---
layout: post
title:  "Hartree-Fock"
date:   2017-05-08
categories: blog
use_math: true
---

When beginning graduate school, I joined Dr. Madura's computational research group with absolutely no programming experience. In implementing computational chemistry methods, density functional theory (DFT) was the first. However, understanding the theory of DFT required me to step back to other methods. Hartree-Fock tends to be the foundation for many other computational chemistry methods.  Reading through Szabo and Ostlund's Modern Quantum Chemistry was helpful, but not quite enough to clarify the algorithm. So what you will find below is some framework code and theory I used to calculate HF. The full code can be found on my [github.](https://github.com/amandadumi/HartreeFock "Amanda's Hartree-Fock repository")

Keeping in mind that our overall goal is to solve the  Schr&ouml;dinger equation,
$$H\Psi = E \Psi$$ we need to establish the form of our wavefunction. Here we will used a single determinant, in other words, a set of spin orbitals that when placed in a determinant represent the ground state energy. The restraint on these orbitals is that they be orthonormal, or in other words
For my code I utilized [PySCF](https://github.com/sunqm/pyscf "PySCF repository") to extract the atomic orbitals that make up my determinant.

$$\langle\Psi\rvert = \chi_1\chi_2 \dots \chi_a\chi_b \dots \chi_N \quad \quad  \langle\chi_a\rvert\lvert\chi_b\rangle =\ \delta_{ab}$$

 The Hamiltonian for Hartree-Fock is as follows:

$$h(1)\ +\ \sum\limits_{b \neq a}\ J_b(1) + \sum\limits_{b \neq a}\ K_b(1)$$

 where $$h(1)$$ will describe the kinetic energy and the potential resulting from the nuclei and an electron:  and $$K_b(1)$$ represents the exchange operator.

  <!-- We can then express the problem as an eigenvalue problem as follows

$$[h(1)\ +\ \sum\limits_{b \neq a}\ J_b(1) + K_b(1)]\chi_a\ = \ \epsilon\chi_a$$ -->

 $$h(1)\ =\ \frac{1}{2} \nabla^2 - \sum\limits_a\frac{Z_a}{r_{12}}$$

$$J_b(1)  =\ \int dx_2 \lvert\chi_b(2)\rvert ^2 r_{12}^{-1}$$ represents the Coulomb operator and describes the force on electron 1 in $$\chi_a$$ affected by electron 2  in $$\chi_b$$. However, we replace the explicit expression by a one electron potential obtained by averaging the interaction $r_{12}^{-1}$. it is said to describe the average local potential at $$x_1$$ arrived from $$\chi_b$$.


 $$K_b(1)\ =\ \int dx_2\chi_b^*(2)r_{12}^{-1}\chi_a(2) $$


First a molecule is defined using PySCF. We will be treating water as an example. We will describe our basis set with sto-3g. The code block below describes the position of each atom in cartesian coordinates.


    mol = gto.M(
    atom = [['O', (0.000000000000,  -0.143225816552,   0.000000000000)],
            ['H', (1.638036840407,   1.136548822547,  -0.000000000000)],  
            ['H', (-1.638036840407,   1.136548822547,  -0.000000000000)]],  
    #basis = '6-31G',  
    basis = 'sto-3g',  
    verbose = 1,  
    unit='b',  
    symmetry=False)  


We can then extract some of the necessary values  



    mo_coefficients = myhf.mo_coeff



#Steps to Hartree-Fock

##### 1. Define the variable for the nuclear-nuclear repulsion energy  


    Enuc = gto.mole.energy_nuc(mol)
    print Enuc


#####2. Read in the one-electron integrals.  

    ovlp = myhf.get_ovlp()
    h1e = myhf.get_hcore()


#####3. Read in the two-electron integrals

    eri = mol.intor('cint2e_sph',aosym=4)
    print np.shape(eri)


#####4. As a check described by Evangelista, check that the following sums are correct

#####5. Form the inverse of the overlap matrix.

        def orthogonalize(Smatrix):
            S_eval,S_evec = spla.eigh(Smatrix)
            #diagonalize the matrix
            S_eval_square=sp.zeros((7, 7),dtype = np.float64)
            for i in range(len(S_eval)):
               S_eval_square[i,i] = S_eval[i]

            eval_change=sp.zeros((len(S_eval),len(S_eval)),dtype = np.float64)
            for i in range(len(eval_change)):
                eval_change[i,i]=1/sp.sqrt(S_eval[i])

            SOM = np.dot(S_evec,np.dot(eval_change,S_evec.T))
            return SOM


#####6. form the G matrix

    def G_matrix(old_Dmatrix,two_e_integrals):
        G = np.zeros(np.shape(old_Dmatrix))
        c=0
        for i in range(len(old_Dmatrix)):
            for j in range(len(old_Dmatrix)):
                for m in range(len(old_Dmatrix)):
                    for n in range(len(old_Dmatrix)):
                        c += (old_Dmatrix[m,n])*((2.0*two_e_integrals[idx4(i,j,m,n)])-two_e_integrals[idx4(i,m,j,n)])
                G[i,j] = c
                c = 0
        return G



#####7. Form the Fock matrix and compute the total energy


    def F_matrix(Hcore,Gmatrix):
            F = np.zeros(np.shape(Hcore))
            for i in range(len(F)):
                for j in range(len(F)):
                    F[i,j] = Hcore[i,j] + G[i,j]
            return F



#####8. now transform the Fock matrix to represent orthonormal atomic orbitals.

    def orthogonalize(Smatrix):
        S_eval,S_evec = spla.eigh(Smatrix)
        #diagonalize the matrix
        S_eval_square=sp.zeros((7,7),dtype = np.float64)
        for i in range(len(S_eval)):
           S_eval_square[i,i] = S_eval[i]

        eval_change=sp.zeros((len(S_eval),len(S_eval)),dtype = np.float64)
        for i in range(len(eval_change)):
            eval_change[i,i]=1/sp.sqrt(S_eval[i])

        SOM = np.dot(S_evec,np.dot(eval_change,S_evec.T))   #This is the symmetric orthogonalization matrix. Right compared to crawdad
        #Y = np.dot(np.transpose(SOM),np.dot(S,SOM)) #should be an identity matrix if done correctly.
        #print Y
        return SOM




#####9. Find the eigenvalue and eigenvectors of the Fock matrix in the orthonormal basis set and transform them back into the atomic orbital basis set.

    def diagonalize(ortho, Fmatrix):
        F_prime = sp.dot(SOM,sp.dot(F,SOM))
        f_eval, f_evec = spla.eigh(F_prime)
        C = sp.dot(SOM,f_evec)
        return C

#####10. Form the density matrix

    def D_matrix(Hcore,C):
        D = np.zeros(np.shape(Hcore))
        temp = 0
        for u in range(len(D)):
            for v in range(len(D)):
                for i in range(num_elec/2):
                    temp += C[u,i] * C[v,i]
                D[u,v]=temp
                temp = 0
        return D



#####11. Test for convergence
  a. Energy Difference
  b. root mean square error of the density.

def RMSD(new_Dmatrix,old_Dmatrix):
        RMS = new_Dmatrix - old_Dmatrix
        RMS = np.power(RMS,2)
        RMS = np.sum(RMS)
        RMS = np.sqrt(RMS)
        return RMS




#####12. Initialize your variables

    myhf = scf.RHF(mol)
    E_hf = myhf.kernel()
    D=np.zeros((len(Hcore),len(Hcore)),dtype = np.float64)
    k=0
    E_threshold = 1e-9
    D_threshold = 1e-5
    convergence = False
    Etot_old  = 0
    Etot_new  = 0




#####14. Repeat until convergence.

    while (convergence == False):
        print '###################### iteration',k,' #######################'
        Dold = np.copy(D)

        ### G matrix ###
        G=G_matrix(Dold,eri)
        k+=1

        ### Fock matrix ###
        F = F_matrix(Hcore,G)

        ###calculate total energy###
        Etot_old = Etot_new
        Eelec = Eelec_calc(Dold,Hcore,F)
        Etot_new = Enuc + Eelec

        #### Find orbital coefficients ###
        C = diagonalize(SOM,F)

        ### Denisty matrix ###
        D = D_matrix(Hcore,C)

        #### RMSD of density matrix ###
        RMS = 0
        RMS=RMSD(D,Dold)
        print 'RMS \t',RMS

        ### Energy difference ###
        E_diff = np.abs(Etot_new - Etot_old)

        if (E_diff < E_threshold) and (RMS < D_threshold):
            convergence = True

        print 'E_diff: ', E_diff

    print 'final (Ha)': {}.format( tot_new)
    print 'error (Ha): {}'.format(np.abs(Etot_new - true_value))
    print 'final step count: {}'.format(k)
