//! Simulated sampling from Markov networks.
//!
//! Simulation of Markov networks probabilities can be difficult because computation of the
//! partition function is often infeasible. Here are some alternative methods that try to
//! simulate by avoiding the partition function.

use ndarray::{Array, Array1, Array2};
use polars::frame::DataFrame;

/// Simulate using a SRBM(sparse restricted Boltzmann model)-based Gibbs sampler.
///
/// # Arguments
/// * `data` - A string slice that holds the name of the person
/// * `coupling_parameters` - An array containing the
/// * `bias_parameters`
pub fn simulate_srbm(
    coupling_parameters: Array2<f64>,
    bias_parameters: Array1<f64>,
    n: usize,
    t: usize,
    b: usize,
) -> DataFrame {
    let mut simulation = DataFrame::from(&data.schema());

    // Introduce hidden variables by creating a new matrix with coupling parameters
    // between the original visible and hidden variables.
    let mut boltzmann_coupling = Array2::from_elem((3, 3), 0);

    // Burn-in
    let mut b = b;
    while b > 0 {

        b -= 1;
    }

    // Simulate
    let mut i = 0;
    let mut j = 0;
    while i < n {

        if j < t {
            j += 1;
        } else {

        };
    };
    Array2::from_diag

    simulation
}
