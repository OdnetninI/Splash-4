Splash-4 Benchmark Suite
========================

Splash-4 is a benchmark suite based on Splash-3 but with enfasis in atomic instructions. It was first presented at ISPASS'21. Then, fully released at IISWC'22.

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

    Eduardo José Gómez-Hernández, Juan Manuel Cebrian, 
    Stefanos Kaxiras, Alberto Ros, "Splash-4: A Modern Benchmark
    Suite with Lock-Free Constructs". 2022 IEEE International
    Symposium on Workload Characterization (IISWC), pages 56--64, 
    Austin, Texas, November 2022

### BibTeX
    
    @InProceedings{ejgomez-iiswc22,
      author =       {Eduardo Jos{\'e} G{\'o}mez-Hern{\'a}ndez and Juan Manuel Cebrian and Stefanos Kaxiras and Alberto Ros},
      title =        {Splash-4: A Modern Benchmark Suite with Lock-Free Constructs},
      booktitle =    {2022 IEEE International Symposium on Workload Characterization (IISWC)},
      doi =          {10.1109/IISWC55918.2022.00015},
      pages =        {51-64},
      year =         {2022},
      editor =       {},
      address =      {Austin, TX (USA)},
      month =        nov,
      publisher =    {IEEE Computer Society},
      ratio-acep =   {47.92% (23/48)},
      isbn =         {},
      url =          {http://webs.um.es/aros/papers/pdfs/ejgomez-iiswc22.pdf}
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
