# RandomizedSparsification.jl

A Julia package implementing the randomized element-wise matrix sparsification algorithm from Kundu and Drineas (2014).

The main function `sparsify(A; r)` sparsifies a matrix `A` by sampling `r` elements according to a probability distribution based on the matrix's Frobenius and 1-norms.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kunzaatko.github.io/RandomizedSparsification.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kunzaatko.github.io/RandomizedSparsification.jl/dev/)
[![Build Status](https://github.com/kunzaatko/RandomizedSparsification.jl/actions/workflows/CI.yml/badge.svg?branch=trunk)](https://github.com/kunzaatko/RandomizedSparsification.jl/actions/workflows/CI.yml?query=branch%3Atrunk)
[![Coverage](https://coveralls.io/repos/github/kunzaatko/RandomizedSparsification.jl/badge.svg?branch=trunk)](https://coveralls.io/github/kunzaatko/RandomizedSparsification.jl?branch=trunk)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
<a href="https://kunzaatko.github.io/">
    <img src="https://raw.githubusercontent.com/pedromxavier/flag-badges/main/badges/CZ.svg" alt="made in Czechia">
</a>

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
