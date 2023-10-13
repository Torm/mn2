//! Discovery of multinomial Markov networks
//!
//! Goal: generate a ranking of arcs.

use std::ops::Index;
use statrs::function::gamma::ln_gamma;
use statrs::statistics::Statistics;

pub mod count;
pub mod iter;
//pub mod format;
//pub mod simulate;
pub mod discovery;

/// Heap array.
pub type Array<T> = Box<[T]>;

/// Represents a sequence of multinomial variables. An index represents a variable
/// number, while a value represents the number of values in the variable range.
pub struct Variables {
    variables: Array<u32>,
}

impl Variables {

    pub fn new(variables: Box<[u32]>) -> Self {
        Self { variables }
    }

    pub fn len(&self) -> usize {
        self.variables.len()
    }

}

impl Index<u32> for Variables {

    type Output = u32;

    fn index(&self, index: u32) -> &Self::Output {
        &self.variables[index as usize]
    }

}

/// Represents a configuration: a variable taking on a specific value.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct Config {
    pub variable: u32,
    pub value: u32,
}

impl Config {

    pub fn new(variable: u32, value: u32) -> Self {
        Self { variable, value }
    }

    pub fn variable(&self) -> u32 {
        self.variable
    }

    pub fn value(&self) -> u32 {
        self.value
    }

}

//// Utility math functions

/// The log-score of a particular variable-parent config depends on the count of all the variables
/// and the count of the parent config.
///
/// - *variable_count* Count of variable and parent config.
/// - *parent_count* Count of parent config.
/// - *variable_alpha* Hyperparameter depending on parent config. `$\alpha_{x_i|u_i}$`
/// - *config_aplha* Hyperparameter depending on variable config and parent config.
///
/// ## K2 prior
///
/// The K2 prior uses a fixed `alpha` value for each child variable.
///
/// Very simple, but is inconsistent: prior observations grows with the number of
/// parents.
///
/// 1) Choose an `alpha` that represents the number of prior observations.
/// 2) `variable_alpha` is equal to `alpha`.
/// 3) `config_alpha` is `alpha` divided by the number of configs of `x`.
///
/// ## BDe (uniform) prior
///
/// The BDe prior uses an `variable_alpha` value for a child variable that depends
/// on the prior probability of each parent config. For a uniform prior over the parent
/// configs, `variable_alpha` can be found by dividing `alpha` by the cardinality
/// of the parent configs.
///
/// 1) Choose an `alpha` that represents the number of prior observations.
/// 2) `variable_alpha` is `alpha` divided by the total number of parent configurations.
/// 3) `config_alpha` is `variable_alpha` divided by number of configs of `x`.
///
pub fn logscore(variable_count: u32, parent_count: u32, variable_alpha: f64, config_alpha: f64) -> f64 {
    ln_gamma(variable_alpha) - ln_gamma(variable_alpha + parent_count as f64) + ln_gamma(config_alpha + variable_count as f64) - ln_gamma(config_alpha)
}

/// See [logscore].
///
/// Logscore when the count of all configs are zero, but the count of the parent config is nonzero.
pub fn logscore_prior(parent_count: u32, variable_alpha: f64) -> f64 {
    ln_gamma(variable_alpha) - ln_gamma(variable_alpha + parent_count as f64)
}

pub fn logscore_config(variable_count: u32, config_alpha: f64) -> f64 {
    ln_gamma(config_alpha + variable_count as f64) - ln_gamma(config_alpha)
}

#[test]
fn ffff() {
    let alpha_v = 1.0/64.0; let alpha_c = 1.0/128.0; let pc = 1.0;
    let prior = ln_gamma(alpha_v) - ln_gamma(alpha_v + pc);
    let c1 = ln_gamma(alpha_c) - ln_gamma(alpha_c);
    println!("Prior[{}] c1[{}] av[{}] ac[{}] r={}", prior, c1, ln_gamma(alpha_v), ln_gamma(alpha_c), prior + c1);
    println!("Loggamma {} = {}", 0.0, ln_gamma(0.0));
    println!("Loggamma {} = {}", 0.1, ln_gamma(0.1));
    println!("Loggamma {} = {}", 0.2, ln_gamma(0.2));
    println!("Loggamma {} = {}", 0.3, ln_gamma(0.3));
    println!("Loggamma {} = {}", 0.5, ln_gamma(0.5));
    println!("Loggamma {} = {}", 0.8, ln_gamma(0.8));
    println!("Loggamma {} = {}", 1.0, ln_gamma(1.0));
    println!("Loggamma {} = {}", 1.5, ln_gamma(1.5));
    println!("Loggamma {} = {}", 2.0, ln_gamma(2.0));
    println!("Loggamma {} = {}", 4.0, ln_gamma(4.0));


}

/// Compute the log-sum-exp of a set of logged values.
///
/// https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
pub fn log_sum_exp(logs: &[f64]) -> f64 {
    if logs.is_empty() {
        return 0.0;
    }
    let max = logs.max();
    log_sum_exp_with_max(logs, max)
}

/// See [log_sum_exp].
pub fn log_sum_exp_with_max(logs: &[f64], max: f64) -> f64 {
    logs.iter()
        .map(|l| (l - max).exp())
        .sum::<f64>()
        .ln()
        + max
}

pub fn log_sum_exp_two(logs_a: &[f64], logs_b: &[f64]) -> f64 {
    if logs_a.is_empty() {
        return log_sum_exp(logs_b);
    };
    if logs_b.is_empty() {
        return log_sum_exp(logs_a);
    }
    let max_a = logs_a.max();
    let max_b = logs_b.max();
    let max = if max_a > max_b {
        max_a
    } else {
        max_b
    };
    let mut sum = 0.0;
    for a in logs_a {
        sum += (a - max).exp();
    }
    for b in logs_b {
        sum += (b - max).exp();
    }
    sum.ln()
}

#[test]
fn test_sum_exp_log() {
    assert_eq!(log_sum_exp(&[1000.0, 1000.0, 1000.0]),        1001.0986122886682);
    assert_eq!(log_sum_exp(&[1000.0, 1000.0, 1000.0, 990.0]), 1001.0986274218635);
}
