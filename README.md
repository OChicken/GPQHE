# GPQHE

A C library doing fully homomorphic encryption under the license of LGPL.

## Architecture

<table>
<tr>
  <td colspan="6" style="text-align: center">enc/dec</td>
</tr>
<tr>
  <td colspan="6" style="text-align: center">RLWE CPAKEM</td>
</tr>
<tr>
  <td colspan="3" style="text-align: center">sample</td>
  <td colspan="3" style="text-align: center">polynomial arithmetics</td>
</tr>
<tr>
  <td style="text-align: center">zero center</td>
  <td style="text-align: center">hwt</td>
  <td style="text-align: center">discrete Gaussian</td>
  <td colspan="3" style="text-align: center">poly-mpi &amp; poly-rns</td>
</tr>
<tr>
  <td colspan="3" style="text-align: center">rng</td>
  <td colspan="2" style="text-align: center">ntt</td>
  <td style="text-align: center">rns</td>
</tr>
<tr>
  <td colspan="3" style="text-align: center"></td>
  <td style="text-align: center">Montgomery</td>
  <td style="text-align: center">Barrett</td>
  <td style="text-align: center">thread</td>
</tr>
</table>

<table>
<tr class="header">
  <th style="text-align: left;"></th>
  <th style="text-align: left;">Used for ...</th>
  <th style="text-align: left;">Keys or Functions</th>
</tr>
<tr>
  <td rowspan="3" style="text-align: left;">Advanced</td>
  <td style="text-align: left;">linear transformation</td>
  <td style="text-align: left;">
    <span class="math inline">gemv(<strong>A</strong>,ct,rk)</span>,
    <span class="math inline">sum(ct,rk)</span>,
    <span class="math inline">idx(ct,<em>i</em>,rk)</span>,
    <span class="math inline">nrm22(ct,rlk,rk,ck)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">nonlinear functions</td>
  <td style="text-align: left;">
    <span class="math inline">exp(ct,rlk)</span>,
    <span class="math inline">ln(ct,rlk)</span>,
    <span class="math inline">sigmoid(ct,rlk)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">Comparison</td>
  <td style="text-align: left;">
    <span class="math inline">inv(ct,rlk,<em>d</em>)</span>,
    <span class="math inline">sqrt(ct,rlk,<em>d</em>)</span>,
    <span class="math inline">cmp(ct<sub>1</sub>,ct<sub>2</sub>,rlk,<em>d</em>,<em>α</em>)</span>
  </td>
</tr>
<tr>
  <td rowspan="5" style="text-align: left;">Primitives</td>
  <td rowspan="2" style="text-align: left;">addition</td>
  <td style="text-align: left;">
    <span class="math inline">add<sub>pt</sub>(ct,pt)</span>,
    <span class="math inline">add(ct<sub>1</sub>,ct<sub>2</sub>)</span>,
    <span class="math inline">neg(ct)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">
    <span class="math inline">sub<sub>pt</sub>(ct,pt)</span>,
    <span class="math inline">sub(ct<sub>1</sub>,ct<sub>2</sub>)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">level switch</td>
  <td style="text-align: left;">
    <span class="math inline">rs<sub>ℓ → ℓ′</sub>(ct)</span>,
    <span class="math inline">moddown<sub>ℓ → ℓ′</sub>(ct)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">multiplication</td>
  <td style="text-align: left;">
    <span class="math inline">mult<sub>pt</sub>(ct,pt)</span>,
    <span class="math inline">mult(ct<sub>1</sub>,ct<sub>2</sub>,rlk)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">automorphism</td>
  <td style="text-align: left;">
    <span class="math inline">rot(ct,<em>r</em>,rk)</span>,
    <span class="math inline">conj(ct,ck)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">Encryptor</td>
  <td style="text-align: left;">encrypt/decrypt</td>
  <td style="text-align: left;">
    <span class="math inline">ct = enc<sub>pk/sk</sub>(pt,pk/sk)</span>,
    <span class="math inline">pt = dec(ct,sk)</span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">Encoder</td>
  <td style="text-align: left;">encode/decode</td>
  <td style="text-align: left;">
    <span class="math inline">pt = ecd(<strong>z</strong>,<em>Δ</em>)</span>,
    <span class="math inline"><strong>z</strong> = dcd(pt,<em>Δ</em>)</span>
  </td>
</tr>
<tr>
<td rowspan="4" style="text-align: left;">KEM</td>
<td style="text-align: left;">public key, secret key</td>
  <td style="text-align: left;">
    <span class="math inline">pk, sk = keypair(ctx)</span>
  </td>
</tr>
<tr>
<td rowspan="3" style="text-align: left;">evaluation keys</td>
<td style="text-align: left;">
  <span class="math inline">rlk = genswk(sk⋅sk,sk)</span>
</td>
</tr>
<tr>
  <td style="text-align: left;">
    <span class="math inline">rk<sub><em>r</em></sub> = genswk(<em>τ</em><sub>5<sup><em>r</em></sup>mod2<em>n</em></sub>(sk),sk), 0 ≤ <em>r</em> &lt; ⋯<em>n</em><sub>slots</sub></span>
  </td>
</tr>
<tr>
  <td style="text-align: left;">
    <span class="math inline">ck = genswk(<em>τ</em><sub>−1</sub>(sk),sk)</span>
  </td>
</tr>
</table>

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
