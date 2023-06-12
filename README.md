# GPQHE

[GPQHE](https://github.com/OChicken/GPQHE) is a reference implementation of CKKS homomorphic encryption scheme. It is a free software under the terms of LGNU v2.1, and its architecture is designed analog to CKKS’s [HEAAN](https://github.com/snucrypto/HEAAN) and two NIST PQC submissions [NewHope](https://newhopecrypto.org/) & [Crystal-Kyber](https://pq-crystals.org/kyber/).

## Description

The virtues of [GPQHE](https://github.com/OChicken/GPQHE) includes

- Except using [Libgcrypt](https://gnupg.org/software/libgcrypt/index.html) (the cryptographic library in Linux kernel) for multi-precision integer (MPI) management, it does not rely on any existing HE libraries;
- Support NTT and MPI-RNS (residual number system) conversion in polynomial multiplication;
- Support basic HE primitives, linear transformation, digits extraction, nonlinear functions and ciphertext comparison;
- Purely written in C, and utilize linked-list for context parameters management;
- Valgrind is used to guarantee that no memory leakage is possible.

The architecture of a typical RLWE (that utilizing RNS in polynomial arithmetics)

![RLWE architecture](doc/rlwe.png)

The standard HE workflow:

![HE workflow](doc/he-workflow.png)

The supported algorithms in GPQHE:

![Supported algorithms in GPQHE](doc/gpqhe.png)

## Package dependency

libgcrypt 1.10. This library is used in Linux kernel, and we utilize its mpi (multi-precision integer) module in GPQHE. If your system does not have it, do `sudo apt install libgcrypt11-dev`.

## How to use

```sh
# step 1: get GPQHE and build
git clone https://github.com/OChicken/GPQHE.git
cd GPQHE
git submodule init
git submodule update
mkdir -p lib
make

# step 2: run tests
cd tests
LD_LIBRARY_PATH=$PWD/../lib:$LD_LIBRARY_PATH
make test-gpqhe
./test-gpqhe enc sk
./test-gpqhe mul pk
```

## Remarks

1. For KAT (known answer tests), we use compiler option `-DSUPERCOP` in `src/Makefile` to deterministic generate random numbers. To fully support true rng, change it to `-DRANDOM`.

2. `tests/polymul.gp` is a [PARI/GP](https://pari.math.u-bordeaux.fr/) script file to verify correctness of polynomial multiplication.

## Acknowledgement

This library is developed upon [HEAAN](https://github.com/snucrypto/HEAAN), [NewHope](https://newhopecrypto.org/) and [Kyber](https://pq-crystals.org/kyber/). We show gratefulness to the developers.
