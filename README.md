# Planet Migration with Nested Sampling

## Overview

This project implements planet migration simulations using nested sampling via [Dynesty](https://github.com/joshspeagle/dynesty). It explores the parameter space efficiently to model planetary migration scenarios and fits observational data using Bayesian inference.

## Features

- Simulates planet migration dynamics (e.g., `swift_symba`)
- Uses Dynesty for dynamic nested sampling to infer posterior distributions
- Supports optional parallelization for performance gains
- Generates corner plots of posterior distributions and sampling diagnostics

## Dependencies

Install the necessary Python packages via pip:

```bash
pip install numpy scipy matplotlib emcee corner dynesty pathos dill
```

You'll also need an installed Fortran compiler (`gfortran`) for the migration integrators (`swift_symba`, etc.).

## File Structure

- `compute_mig_nest.py`: Main script for running migration nested sampling  
- `dynesty_2_0/`: Forked Dynesty sampler code  
- `corner_ES/`: Forked cornerplot code  
- `swift/`: Compiled Fortran integrator libraries (e.g. `libswift.a`)  
- `README.md`: This documentation

## Usage

Make sure the swift binaries are compiled

```bash
bash install_swift.sh
```

Run the main sampling script:

```bash
python3 compute_mig_nest.py
```

By default, this will execute migration models with nested sampling. Inside compute_mig_nest.py you can configure:

- Input parameters via command-line flags or configuration files  
- Sampling options (number of live points, stopping criteria)  
- Parallel settings using `pathos.ProcessPool`

Output includes posterior sample files and optional corner plots.

## References

- Speagle, J. S. (2020). *DYNESTY: a dynamic nested sampling package for estimating Bayesian posteriors and evidences*. MNRAS, 493(3), 3132–3158. doi:10.1093/mnras/staa278  
- [Swift integrator code documentation](http://www.astro.umd.edu/~hamilton/research/swift.html)

## License

Open-source under the MIT License. Contributions and issues are welcome.

## About

Developed by Trifon Trifonov. 
