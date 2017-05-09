---
layout: post
title:  "Hartree-Fock"
date:   2017-05-08
---

<!-- When beginning graduate school, I joined Dr. Madura's computational research group with absolutely no programming experience. In learning about computational chemistry methods, one of the first methods I learned about and implemented was density functional Theory. However, my understanding of DFT and other computational chemistry methods wasn't really in perspective until  I learned about Hartree-Fock. Reading through Szabo and Ostlund was helpful, but not quite enough to clarify the algorithm. So what you will find below is some framework code and theory I used to calculate HF. -->

A full copy of the code can be found on my github through this [link](https://github.com/amandadumi/HartreeFock)

Keeping in mind that our overall goal is to solve the  Schr&ouml;dinger equation, $H\Psi = E \Psi$, we need to establish the form of our wavefunction. Here we will used a single determinant, or in other words a set of spin orbitals, when placed in a determinant represent the ground state energy.
$$\Psi_0$$
