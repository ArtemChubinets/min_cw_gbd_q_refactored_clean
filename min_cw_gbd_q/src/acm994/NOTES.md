# ACM CALGO 994 minimum-codeword adapter

This adapts [the official ACM CALGO 994 release](https://calgo.acm.org/994.zip),
*Fast Implementations of the Brouwer-Zimmermann Algorithm for the Computation
of the Minimum Distance of a Random Linear Code*, Fernando Hernando,
Francisco D. Igual, and Gregorio Quintana-Ortí, ACM TOMS 45(2), Article 23,
June 2019. The archive SHA-256 is
`f08357bcbd9b48abd42107a7619b3e07d95986a492f7ce7a77567f03e4c46d72`.

The upstream source headers explicitly grant GPL-3.0-or-later; the archive's
`994/Doc/LICENSE` contains GPL v3. The build preserves the upstream notices
and installs that license under `share/licenses/acm994-codeword/LICENSE`.
The adapter is also GPL-3.0-or-later. The [CALGO index](https://calgo.acm.org/)
points to the [ACM software notice](https://www.acm.org/publications/policies/softwarecrnotice).
This repository keeps the adapter and patch, not the upstream archive,
upstream source tree, or compiled executable. If distributing the executable,
ensure corresponding source and applicable notices are provided as required;
prefer the local build below.

## Build and run

On Linux/WSL, with Python 3, `curl`, GNU `patch` and an OpenMP-capable C
compiler, run from the repository root:

```sh
CC=gcc python3 min_cw_gbd_q/src/acm994/build.py
```

Add `--archive /path/to/994.zip` for an offline build, or `--prefix /path/to/install`
to change the installation directory. Both fetch modes verify the same hash.
For a host-local optimized build, set `CFLAGS='-O3 -march=native'` explicitly.
The default is portable scalar `-O3` with GCC; `-march=native` can use POPCNT
on a supporting build host. The compiler in the local WSL Sage environment
can be invoked as follows:

```sh
/home/artem/.local/share/qgbd-tools/bin/micromamba run \
  -r /home/artem/.local/share/qgbd-tools/mamba -n sage \
  env CFLAGS='-O3 -march=native' python min_cw_gbd_q/src/acm994/build.py
```

The default executable is `~/.local/share/qgbd-tools/acm994/bin/acm994-codeword`.
It accepts `k n` followed by `k` lines of `n` contiguous binary digits, and
prints one stdout line `RESULT <weight> <n-bit word>`. The word is in input
column order. Dependent rows are reduced first; a zero code or malformed
input exits nonzero without `RESULT`. Set `ACM994_VERBOSE=1` for search logs
on stderr.

```sh
printf '2 4\n1010\n1100\n' | ~/.local/share/qgbd-tools/acm994/bin/acm994-codeword
```

## What changes from upstream

- Retains CALGO 994's Gamma diagonalization, permutation selection (five
  deterministic trials) and saved-combination search. Uses the serial scalar
  `ALG_SAVED` variant, caching sizes one and two (one only at rank one).
- Full packed rows include pivot coordinates, rather than dropping the
  identity portion and adding its weight back. The saved kernel keeps the
  prefix and winning suffix for its minimum; the hook XORs them and reverses
  the selected Gamma's column order if they improve the global upper bound.
  There is no second search or code shortening. The CLI also checks the
  resulting word's weight and membership in the original row span.
- The loop visits every Gamma at each coefficient weight until the bounds
  meet; it does not use upstream's contributing-Gamma pruning. For a Gamma
  of rank `r` searched through coefficient weight `i`, a word not yet found
  uses at least `i+1` rows, so at least `max(0, i+1-(k-r))` of its pivot bits
  are ones, where `k` is the input rank. The Gaussian column swaps are confined
  to the trailing columns, so the actual processed pivot blocks are disjoint
  and these lower bounds can be added. A retained codeword is an upper bound;
  equality or crossing proves its minimality. If weight `k` is searched in
  any basis, that basis alone has enumerated every nonzero word.
- These choices change the cost from the unmodified upstream variants:
  more coordinates can be counted per candidate, more Gamma matrices can
  be visited, and the search does not use multiple workers. Timings of this
  executable should not be labeled as upstream parallel/SIMD or posted
  CALGO timings.
