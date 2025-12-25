**GameTalker**
GameTalker is a computational framework for modeling and quantifying
dynamic intercellular interactions in multicellular systems using
time-varying ordinary differential equation (ODE) models with
functional interaction terms.

The framework is designed to disentangle intrinsic growth dynamics and
interaction-driven effects, and to infer how pairwise and higher-order
interactions evolve over time from longitudinal abundance or expression data.

This repository provides the full implementation of the GameTalker model,
together with a minimal reproducible example illustrating the complete
analysis workflow used in the accompanying manuscript.

## Repository Structure

```text
GameTalker/
├── main.R                 # Core model implementation
├── example_inference.R    # Minimal reproducible example
├── README.md
```

File descriptions

main.R
Implements the complete GameTalker framework, including:
intrinsic growth models,
pairwise and higher-order interaction ODEs,
time-varying interaction functions,
numerical solvers (Runge–Kutta),
objective functions and parameter inference utilities.

example_inference.R
A self-contained toy example demonstrating the standard GameTalker workflow.
