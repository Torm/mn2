//! Discovery of Markov network structures

use std::slice::Iter;

/// A ranking.
///
/// This is the data structure returned by the ranking function. It is a list of
/// arcs ordered by which arc increases the score of the Markov network most on average,
/// within the search space.
pub struct Ranking {
    entries: Vec<(u32, u32, f64)>,
}

impl Ranking {

    /// Create a ranking from a sorted vector.
    ///
    /// **Does not check that the vector is sorted!**
    pub fn from_sorted_vec(entries: Vec<(u32, u32, f64)>) -> Self {
        Ranking { entries }
    }

    pub fn iter(&self) -> Iter<(u32, u32, f64)> {
        self.entries.iter()
    }

}
// todo: Main entry function
// //// Rank arcs.
// //// Main method with many configuration options.
//
// /// Create a ranking of arcs given the data.
// ///
// ///
// ///
// /// Given `m` observations of these variables, create a ranking of which arcs are most likely
// /// included.
// ///
// /// We do this by computing
// ///
// pub fn rank_arcs(input: InputMethod, method: RankMethod, merge_strategy: MergeMethod) -> Result<Ranking, String> {
//     let log_scores = match method {
//         RankMethod::KNeighbours(k) => todo!(),
//         RankMethod::MonteCarlo => todo!(),
//     };
// }

//// Ranking method

/// Ranking method.
pub enum RankMethod {
    KNeighbours(u32),
    MonteCarlo,
}

//// Score merging

/// Which method to use for merging of scores.
///
/// For example, the affinity `A` has for `B` might be different from the affinity
/// `B` has for `A`. The score merging method determines how to merge these scores
/// to produce a final score for the affinity between `A` and `B`.
pub enum MergeMethod {
    Average,
}

fn merge_scores(method: MergeMethod, a: f64, b: f64) -> f64 {
    match method {
        MergeMethod::Average => (a + b) / 2.0,
    }
}

//// Log score accumulation

/// Criterion for accumulation of logscores.
///
/// Comes down to a tradeoff between memory and accuracy. However, a majority of log-scores are insignificant
/// in comparison to the greatest log-score, thus, there is an opportunity to save massive amounts of
/// space at zero, or almost zero, cost.
pub enum LogScoreAccumulationMethod {
    /// Accumulate all logscores.
    All,
    /// Discard all logscores significantly lesser than the maximum logscore.
    ///
    /// This method discards a majority of logscores assumed that only a percentage
    /// of them will be significant enough.
    Threshold,
}

/// k-neighbours based discovery methods
///
/// These methods consider for each variable all sets of up to k neighbours.
///
/// For every variable x, all k-combinations of possible neighbours are considered.
/// For each such neighbour set, the logscore of how likely these variables are to
/// be neighbours of x is computed. That computed logscore is stored in a collection
/// accounting the logscores from `y` to `x`, for each `y` in the neighbour set.
///
/// To compute a final score for the affinity of the arc from `y` to `x`, all those
/// log-scores must be exponentiated and thereafter summed.
///
/// These are the steps needed to compute a score for each arc:
/// 1) Count the configurations in the input data.
/// 2) Compute and accumulate partial family scores.
/// 3) Compute the family scores.
/// 3) Compute asymmetric scores.
/// 4) Merge asymmetric scores, using some merge method.
/// 5) Create a ranking of arcs, sort by score.
pub mod k_neighbours {
    use std::cmp::max;
    use std::fs::File;
    use std::sync::Arc;
    use std::thread;
    use dashmap::DashMap;
    use fxhash::FxBuildHasher;
    use parking_lot::{Mutex, RwLock};
    use rayon::slice::ParallelSliceMut;
    use log::log;
    use crate::{Array, Config, logscore, log_sum_exp, log_sum_exp_two, Variables, logscore_prior, logscore_config};
    use crate::count;
    use crate::discovery::{Ranking};
    use crate::iter::{ChildConfigIter, VariableIter};

    pub enum InputMethod {
        File(String),
    }

    pub enum LogScoreOption {
        Matrix,
    }

    pub enum ThresholdOption {
        All,
        Threshold(f64),
    }

    // todo: Main function with configurations
    // pub fn score_arcs_by_k_neighbours(
    //     data: InputMethod,
    //     k: u32,
    //     count_method: CountOption,
    //     averaging: MergeMethod,
    //     log_score_option: LogScoreOption,
    //     log_score_threshold: ThresholdOption,
    // ) -> Result<Ranking, String> {
    //     match count_method {
    //         CountOption::Hashed { par } => {
    //             let file = match File::open(data) {
    //                 Ok(f) => f,
    //                 Err(e) => return Err(format!("Error opening file: {}", e.to_string())),
    //             };
    //             let (variables, configs) = if par {
    //                 count::count_file_configs_into_map_par(file, k as usize)?
    //             } else {
    //                 count::count_file_configs_into_map(file, k as usize)?
    //             };
    //             let log_scores = match log_score_option {
    //                 LogScoreOption::Matrix => {
    //                     let log_scores = vec![vec![], n * n];
    //                     match log_score_threshold {
    //                         ThresholdOption::All => {
    //                             all_logscores_from_map()
    //                         }
    //                         ThresholdOption::Threshold(logdiff) => {
    //                             let maximum;
    //                             thresholded_logscores_from_map()
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         CountOption::Tree => score_arcs_by_k_neighbours_tree(),
    //     }
    //     todo!()
    // }

    /// *Implements steps: 1*
    pub fn count_map_from_file(file: File, k: u32) -> Result<(Variables, DashMap<Array<Config>, u32, FxBuildHasher>), String> {
        count::count_file_configs_into_map(file, k as usize + 1)
    }

    /// *Implements steps: 1*
    pub fn count_map_from_file_par(file: File, k: u32) -> Result<(Variables, DashMap<Array<Config>, u32, FxBuildHasher>), String> {
        count::count_file_configs_into_map_par(file, k as usize + 1, None)
    }

    /// A logscore matrix holds a [Vec] of logscore values for each of the possible
    /// arcs `y -> x` between the nodes.
    pub struct LogscoreMatrix {
        n: usize,
        /// A logscore will be discarded, rather than appended, if it is lesser than
        /// `max - threshold`.
        threshold: f64,
        arcs: Vec<LogscoreMatrixArc>,
    }

    #[derive(Clone)]
    pub struct LogscoreMatrixArc {
        /// The maximum found family score.
        max: f64,
        /// A list of logscores. The flag is true if y is in the parent set of x, otherwise it is false.
        logscores: Vec<(bool, f64)>,
    }

    impl LogscoreMatrix {

        /// Create a new logscore matrix.
        fn new(n: u32, threshold: Option<f64>) -> Self {
            let n = n as usize;
            Self {
                n,
                threshold: threshold.unwrap_or(f64::MAX),
                arcs: vec![LogscoreMatrixArc { max: f64::MIN, logscores: vec![] }; n * n],
            }
        }

        ///
        ///
        /// If the flag is true, this logscore was computed where p was in the parent set of x. Otherwise,
        /// p was not in the parent set.
        fn append(&mut self, p: u32, x: u32, flag: bool, logscore: f64) -> bool {
            let threshold = self.threshold;
            let arc = self.get_mut(p, x);
            if arc.max - threshold > logscore {
                return false; // Do not append if it is too small compared to the largest found value.
            }
            if logscore > arc.max {
                arc.max = logscore;
            }
            arc.logscores.push((flag, logscore));
            true
        }

        fn get(&self, p: u32, x: u32) -> &LogscoreMatrixArc {
            self.arcs.get(self.n * p as usize + x as usize).unwrap()
        }

        fn get_mut(&mut self, p: u32, x: u32) -> &mut LogscoreMatrixArc {
            self.arcs.get_mut(self.n * p as usize + x as usize).unwrap()
        }

    }

    /// Compute the family scores for every variable and combination of up to k parent variables.
    ///
    /// - Implements steps: 2
    ///
    /// Returns (logscore_matrix, max, min, accumulated, discarded)
    pub fn familyscores_from_counts(k: u32, variables: &Variables, counts: DashMap<Array<Config>, u32, FxBuildHasher>, threshold: Option<f64>) -> Result<(DashMap<(u32, Array<u32>), f64, FxBuildHasher>, u32), String> {
        let family_scores: DashMap<(u32, Array<u32>), f64, FxBuildHasher> = DashMap::default();
        let mut accumulated = 0; // Number of partial family score terms that were accumulated.
        let counts = counts.into_read_only();
        for (parent_config, parent_config_count) in counts.iter() {
            let parent_len = parent_config.len();
            if parent_len == k as usize + 1 {
                continue; // There are at most k parent variables.
            }
            let mut parent_cardinality = 1;
            for p in parent_config.iter() { // Compute the number of parent configs.
                parent_cardinality *= variables[p.variable] + 1;
            }
            let alpha = 0.5;
            let variable_alpha = alpha / parent_cardinality as f64;
            // Iterate over all child variables and all the child variable's configs.
            let mut iter = ChildConfigIter::new(variables, parent_config);
            iter.next();
            loop {
                let current_child = iter.child().variable;
                let mut partial_family_score = logscore_prior(*parent_config_count, variable_alpha);
                let config_alpha = variable_alpha / iter.child_cardinality() as f64;
                loop { // Iterate over the child variable's configs.
                    if let Some(child_count) = counts.get(iter.buffer()) {
                        partial_family_score += logscore_config(*child_count, config_alpha);
                    }
                    iter.next();
                    if iter.child().variable != current_child {
                        break;
                    }
                }
                accumulated += 1;
                let mut parents = Vec::with_capacity(parent_config.len());
                for pc in parent_config.iter() {
                    parents.push(pc.variable);
                }
                // Add the partial family score term to the P -> x family score entry, where P is the
                // parent variables, and x is the current child variable.
                if let Some(mut family_score) = family_scores.get_mut(&(current_child, parents.clone().into_boxed_slice())) {
                    *family_score += partial_family_score;
                } else {
                    family_scores.insert((current_child, parents.clone().into_boxed_slice()), partial_family_score);
                }
                if iter.is_done() {
                    break;
                }
            }
        }
        Ok((family_scores, accumulated))
    }

    /// Process the family scores.
    ///
    /// For each family score, iterate over every variable y. If this variable is in the parent set, this
    /// family score is added to the positive and the total lists of y -> x. If this variable is not in
    /// the parent set, this family score is added only to the total list of y -> x.
    ///
    ///
    pub fn process_family_scores(n: u32, family_scores: DashMap<(u32, Array<u32>), f64, FxBuildHasher>, threshold: Option<f64>) -> Result<(LogscoreMatrix, u64, u64), String> {
        let mut logscores = LogscoreMatrix::new(n, threshold);
        let mut accumulated = 0;
        let mut discarded = 0;
        let family_scores = family_scores.into_read_only();
        for ((child, parents), score) in family_scores.iter() {
            let mut p = 0;
            loop {
                if p == *child { // todo efficiency
                    p += 1;
                    if p >= n {
                        break;
                    }
                    continue;
                }
                let appended = if parents.contains(&p) {
                    logscores.append(p, *child, true, *score)
                } else {
                    logscores.append(p, *child, false, *score)
                };
                if appended {
                    accumulated += 1;
                } else {
                    discarded += 1;
                }
                p += 1;
                if p >= n {
                    break;
                }
            }
        }
        Ok((logscores, accumulated, discarded))
    }

    /// The log-sum-exp method is used here.
    pub fn process_logscore_matrix(n: u32, logscores: LogscoreMatrix) -> Result<Ranking, String> {
        let mut ranking = vec![];
        for i in 0..n {
            for j in i+1..n {
                // Probability of i as a parent of j
                let s1 = logscores.get(i, j);
                let mut s1_positive = 0.0;
                let mut s1_negative = 0.0;
                for (flag, score) in &s1.logscores {
                    if *flag {
                        s1_positive += (score - s1.max).exp();
                    } else {
                        s1_negative += (score - s1.max).exp();
                    }
                }
                let s1_pos = s1_positive.ln() + s1.max;
                let s1_total = (s1_positive + s1_negative).ln() + s1.max;
                let log_p1 = s1_pos - s1_total;
                // Probability of j as a parent of i
                let s2 = logscores.get(j, i);
                let mut s2_positive = 0.0;
                let mut s2_negative = 0.0;
                for (flag, score) in &s2.logscores {
                    if *flag {
                        s2_positive += (score - s2.max).exp();
                    } else {
                        s2_negative += (score - s2.max).exp();
                    }
                }
                let s2_pos = s2_positive.ln() + s2.max;
                let s2_total = (s2_positive + s2_negative).ln() + s2.max;
                let log_p2 = s2_pos - s2_total;
                // Model averaging of the two probabilities. Could be done by max, min, or mean, for example.
                ranking.push((i, j, log_p1.max(log_p2)));
            }
        }
        ranking.par_sort_unstable_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        Ok(Ranking::from_sorted_vec(ranking))
    }

    // /// todo: Test different locks:
    // /// std Mutex
    // /// parking_lot Mutex
    // /// spin Lock
    // pub struct ConcurrentLogscoreMatrix {
    //     n: usize,
    //     threshold: f64,
    //     max: RwLock<f64>,
    //     arcs: Vec<Arc<Mutex<Vec<f64>>>>,
    // }
    //
    // impl ConcurrentLogscoreMatrix {
    //
    //     fn new(n: u32, threshold: Option<f64>) -> Self {
    //         let n = n as usize;
    //         let mut arcs = Vec::with_capacity(n * n);
    //         for _ in 0..n*n {
    //             arcs.push(Arc::new(Mutex::new(vec![])));
    //         }
    //         Self {
    //             n,
    //             threshold: threshold.unwrap_or(f64::MAX),
    //             max: RwLock::new(f64::MIN),
    //             arcs,
    //         }
    //     }
    //
    //     fn append(&self, p: u32, x: u32, logscore: f64) {
    //         let max = *self.max.read();
    //         if max - self.threshold < logscore {
    //             return;
    //         }
    //         if logscore > max {
    //             let mut locked_max = self.max.write();
    //             if logscore > *locked_max {
    //                 *locked_max = logscore;
    //             }
    //         }
    //         let mut locked_logscores = self.arcs[self.n * p as usize + x as usize].lock();
    //         locked_logscores.push(logscore);
    //     }
    //
    // }

    /// todo Iterate in a different order that can be taken advantage of by multiple threads.
    fn ranking_from_logscore_matrix_par2(_logscores: LogscoreMatrix) -> Result<Ranking, String> {
        todo!()
    }

}

/// Markov-chain Monte-Carlo based discovery methods
pub mod mcmc {

    pub fn rank_mcmc() {
        todo!()
    }

}
