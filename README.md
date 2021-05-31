Splash-4 Benchmark Suite
========================

Splash-4 is a benchmark suite based on Splash-3 but with enfasis in atomic instructions. It was first presented at ISPASS'21.

## Building and Running

Splash-4 does not require any special installation step.

In order to build Splash-4, a relatively modern C compiler with C11 support is
required. In addition, the M4 macro processor is necessary for generating the
source files in the first place. On Ubuntu and similar systems, you can install
the `m4` and `ivtools-dev` packages.

Some benchmarks expect their inputs to be at very specific paths relative to the
working directory, so it is recommended to change the working directory to the
benchmark folder before executing it.

## Citing Splash-4

    Eduardo José Gómez-Hernández, Ruixiang Shao, Christos Sakalis, 
    Stefanos Kaxiras, Alberto Ros, "Splash-4: Improving Scalability 
    with Lock-Free Constructs". International Symposium on Performance 
    Analysis of Systems and Software (ISPASS), pages 235--236, 
    Worldwide event, March 2021

### BibTeX
    
    @InProceedings{ejgomez-ispass21,
      author = 	     {Eduardo Jos{\'e} G{\'o}mez-Hern{\'a}ndez and Ruixiang Shao and Christos Sakalis and Stefanos Kaxiras and Alberto Ros},
      title = 	     {Splash-4: Improving Scalability with Lock-Free Constructs},
      booktitle =    {International Symposium on Performance Analysis of Systems and Software (ISPASS)},
      doi =          {10.1109/ISPASS51385.2021.00044},
      editor =       {IEEE Computer Society},
      pages = 	     {235--236},
      year = 	     {2021},
      address =      {Worldwide event},
      month = 	     mar,
      publisher =    {IEEE Computer Society},
      ratio-acep =   {36.92% (24/65)},
      isbn =         {978-1-7281-8643-6},
    }



## Recommended Inputs

The following are the recommended inputs for running Splash-3 in a simulator.
The symbol '#' is a placeholder for the number of threads. These inputs are
based on the original Splash-2 characterization paper by Woo et al. [1].

	APPS:
	./BARNES < inputs/n16384-p#
	./FMM < inputs/input.#.16384
	./OCEAN -p# -n258
	./RADIOSITY -p # -ae 5000 -bf 0.1 -en 0.05 -room -batch
	./RAYTRACE -p# -m64 inputs/car.env
	./VOLREND # inputs/head 8
	./WATER-NSQUARED < inputs/n512-p#
	./WATER-SPATIAL < inputs/n512-p#
	KERNELS:
	./CHOLESKY -p# < inputs/tk15.O
	./FFT -p# -m16
	./LU -p# -n512
	./RADIX -p# -n1048576

## Known Issues

* If OCEAN (either version) segfaults during initialization, try reducing the
  `IMAX` and `JMAX` values in `decs.h`.
* It is recommended to change the working directory to the benchmark's
  directory.

---

[1] Steven Cameron Woo, Moriyoshi Ohara, Evan Torrie, Jaswinder Pal Singh, and
Anoop Gupta. 1995. The SPLASH-2 programs: characterization and methodological
considerations. In Proceedings of the 22nd annual international symposium on
Computer architecture (ISCA '95). ACM, New York, NY, USA, 24-36.
DOI=http://dx.doi.org/10.1145/223982.223990 
