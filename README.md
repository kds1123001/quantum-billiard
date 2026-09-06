# Quantum Billiards — Terminal

A single-particle 2D quantum mechanics simulation, solved exactly and rendered live in your terminal. No GUI, no browser, no external dependencies beyond LAPACK/BLAS.

## What it actually does

- Builds a finite-difference Hamiltonian on a 2D grid (Dirichlet boundary conditions = infinite walls at the domain edge).
- Obstacles are Gaussian potential bumps placed inside the box — this is what a wavepacket scatters off.
- Diagonalizes the Hamiltonian exactly with LAPACK's `zheev` (one-time cost).
- Projects an initial Gaussian wavepacket onto that eigenbasis.
- Every frame is computed as `psi(t) = sum_n c_n * exp(-i*E_n*t) * phi_n`, evaluated with a single BLAS `zgemv` call. No time-stepping error — evolution is unitary to machine precision (verified: probability norm stays at 1.000000 across a full run).

This is a real quantum mechanics simulation of one particle in a 2D potential landscape — not a many-body electron simulator, not a general-purpose physics engine. That's the honest scope.

## Requirements

sudo apt install -y gfortran liblapack-dev libblas-dev


## Build

make


## Run

./qbilliard_term


You'll be prompted for:

| Prompt | Meaning | Default |
|---|---|---|
| obstacle layout | 0 = empty box, 1 = single bumper, 2 = five bumpers | 2 |
| wavepacket start x/y | initial position of the particle | 4.0 / 5.0 |
| momentum kx/ky | initial momentum (higher = faster) | 6.0 / 0.5 |
| sigma | initial wavepacket width | 0.7 |
| number of frames | how long the animation runs | 400 |
| dt | time step per frame | 0.006 |
| frame delay (ms) | playback speed | 35 |

Press Enter on any prompt to accept the default.

Diagonalization (the one-time setup cost) takes a few seconds on one core — that's genuine `zheev` compute on a 1152-state Hamiltonian, not a fake loading bar. After that, animation is instant per frame.

When a run finishes, it asks if you want to run again with new parameters (`y`/`n`).

## Notes

- Needs a terminal with 24-bit ANSI color and UTF-8 support (Windows Terminal, most Linux terminals, iTerm2 all work fine).
- Grid resolution is set at compile time in `src/geometry.f90` (`nx`, `ny`, `lx`, `ly`). Bigger grid = better resolution but diagonalization time scales as N³, so double the grid points and you're looking at ~8x the setup time.
- Obstacle positions/sizes for each preset are also in `src/geometry.f90`, in `set_obstacles`.

## Project layout

src/kinds.f90 - floating point kind definitions
src/geometry.f90 - grid, obstacle potential, runtime obstacle presets
src/hamiltonian.f90 - Hamiltonian assembly + zheev diagonalization
src/evolve.f90 - wavepacket init + eigenbasis time evolution
src/render.f90 - ANSI terminal rendering + frame timing
src/main.f90 - interactive CLI entry point
Makefile

# credits

all credits go to me yo boi kds1123001 aka quantanamokid46
