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
/// 2) Compute and accumulate logscores.
/// 3) Compute asymmetric scores.
/// 4) Merge asymmetric scores, using some merge method.
/// 5) Create a ranking of arcs, sort by score.
pub mod k_neighbours {

    use std::fs::File;
    use std::sync::Arc;
    use std::thread;
    use dashmap::DashMap;
    use fxhash::FxBuildHasher;
    use parking_lot::{Mutex, RwLock};
    use rayon::slice::ParallelSliceMut;
    use crate::{Array, Config, logscore, log_sum_exp, log_sum_exp_two, Variables, logscore_prior};
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

    /// *Implements steps: 1*
    pub fn score_arcs_by_k_neighbours_tree() -> Result<DashMap<Array<Config>, f64>, String> {
        todo!();
    }

    /// A logscore matrix holds a [Vec] of logscore values for each of the `n * (n - 1)`
    /// arcs `y -> x` between the nodes.
    pub struct LogscoreMatrix {
        n: usize,
        /// A logscore will be discarded, rather than appended, if it is lesser than
        /// `max - threshold`.
        threshold: f64,
        max: f64,
        arcs: Vec<Vec<f64>>,
    }

    impl LogscoreMatrix {

        /// Create a new logscore matrix.
        fn new(n: u32, threshold: Option<f64>) -> Self {
            let n = n as usize;
            Self {
                n,
                threshold: threshold.unwrap_or(f64::MAX),
                max: f64::MIN,
                arcs: vec![vec![]; n * n],
            }
        }

        fn append(&mut self, p: u32, x: u32, logscore: f64) -> bool {
            if self.max - self.threshold > logscore {
                return false;
            }
            if logscore > self.max {
                self.max = logscore;
            }
            let logscores = self.arcs.get_mut(self.n * p as usize + x as usize).unwrap();
            logscores.push(logscore);
            true
        }

    }

    /// - Implements steps: 2
    ///
    /// Returns (logscore_matrix, max, min, accumulated, discarded)
    pub fn logscores_from_map(k: u32, variables: &Variables, counts: DashMap<Array<Config>, u32, FxBuildHasher>, threshold: Option<f64>) -> Result<(LogscoreMatrix, f64, f64, usize, usize), String> {
        let n = variables.len() as u32;
        let mut logscores = LogscoreMatrix::new(n, threshold);
        let mut max = f64::MIN; // The maximum logscore computed.
        let mut min = f64::MAX; // The minimum logscore computed.
        let mut accumulated = 0; // Number of logscores that were accumulated into a LogscoreMatrix.
        let mut discarded = 0; // Number of logscores that were too small to be accumulated.
        let counts = counts.into_read_only();
        for (parent_config, parent_config_count) in counts.iter() {
            let parent_len = parent_config.len();
            if parent_len == k as usize + 1 { // There are at most k parents.
                continue;
            }
            let mut parent_cardinality = 1;
            for p in parent_config.iter() { // Compute the number of parent configs.
                parent_cardinality *= variables[p.variable] + 1;
            }
            let alpha = 0.5;
            let variable_alpha = alpha / parent_cardinality as f64;
            let mut iter = ChildConfigIter::new(variables, parent_config);
            while iter.next() {
                let logscore = if let Some(config_count) = counts.get(iter.buffer()) {
                    let child_cardinality = iter.child_cardinality();
                    let config_alpha = variable_alpha / child_cardinality as f64;
                    logscore(*config_count, *parent_config_count, variable_alpha, config_alpha)
                } else {
                    logscore_prior(*parent_config_count, variable_alpha)
                };
                if logscore < min {
                    min = logscore;
                }
                if logscore > max {
                    max = logscore;
                }
                for Config { variable: parent, .. } in parent_config.iter() { // Add the logscore for u -> x to each p -> x arc.
                    let Config { variable: child, .. } = iter.child();
                    if logscores.append(*parent, child, logscore) {
                        accumulated += 1;
                    } else {
                        discarded += 1;
                    }
                }
            }
        }
        Ok((logscores, max, min, accumulated, discarded))
    }

    /// todo: Test different locks:
    /// std Mutex
    /// parking_lot Mutex
    /// spin Lock
    pub struct ConcurrentLogscoreMatrix {
        n: usize,
        threshold: f64,
        max: RwLock<f64>,
        arcs: Vec<Arc<Mutex<Vec<f64>>>>,
    }

    impl ConcurrentLogscoreMatrix {

        fn new(n: u32, threshold: Option<f64>) -> Self {
            let n = n as usize;
            let mut arcs = Vec::with_capacity(n * n);
            for _ in 0..n*n {
                arcs.push(Arc::new(Mutex::new(vec![])));
            }
            Self {
                n,
                threshold: threshold.unwrap_or(f64::MAX),
                max: RwLock::new(f64::MIN),
                arcs,
            }
        }

        fn append(&self, p: u32, x: u32, logscore: f64) {
            let max = *self.max.read();
            if max - self.threshold < logscore {
                return;
            }
            if logscore > max {
                let mut locked_max = self.max.write();
                if logscore > *locked_max {
                    *locked_max = logscore;
                }
            }
            let mut locked_logscores = self.arcs[self.n * p as usize + x as usize].lock();
            locked_logscores.push(logscore);
        }

    }

    // /// See [logscores_from_map].
    // pub fn logscores_from_map_par(n: u32, counts: DashMap<Array<Config>, u32, FxBuildHasher>, threshold: Option<f64>) -> Result<LogscoreMatrix, String> {
    //     let logscores = ConcurrentLogscoreMatrix::new(n, threshold);
    //     let counts = counts.into_read_only();
    //     //counts.par_iter().for_each(|entry| {
    //     counts.par_iter().for_each(|entry| {
    //         let len = entry.key().len();
    //         if len == 1 { // Skip len 1, since they have no children.
    //             return;
    //         }
    //         let variables = entry.key();
    //         let variables_count = *entry.value();
    //         let mut parents = vec![Config { variable: 0, value: 0 }; len - 1];
    //         for x_index in 0..len {
    //             let x = variables[x_index].variable;
    //             for v in 0..len { // Put all variables, except child, in the parent set
    //                 if v < x_index {
    //                     parents[v] = variables[v];
    //                 } else if v > x_index {
    //                     parents[v - 1] = variables[v];
    //                 }
    //             }
    //             let parents_count = match counts.get((&parents).as_slice()) {
    //                 None => continue,
    //                 Some(l) => *l,
    //             };
    //             let logscore = log_score(variables_count, parents_count);
    //             for p in &parents {
    //                 let p = p.variable;
    //                 logscores.append(p, x, logscore);
    //             }
    //         }
    //     });
    //     let mut arcs = Vec::with_capacity(n as usize * n as usize);
    //     for v in logscores.arcs.into_iter() {
    //         arcs.push(Arc::into_inner(v).unwrap().into_inner());
    //     }
    //     let logscores = LogscoreMatrix {
    //         n: logscores.n,
    //         threshold: logscores.threshold,
    //         max: logscores.max.into_inner(),
    //         arcs,
    //     };
    //     Ok(logscores)
    // }

    /// Compute ranking from logscore matrix.
    ///
    /// *Implements steps: 3, 4, 5*
    pub fn ranking_from_logscore_matrix(logscores: LogscoreMatrix) -> Result<Ranking, String> {
        let n = logscores.n;
        let mut arc_iter = VariableIter::new(n, 2);
        let mut merge = vec![];
        loop {
            if !arc_iter.next() {
                break;
            }
            let i = arc_iter.at()[0] as usize;
            let j = arc_iter.at()[1] as usize;
            let l1 = logscores.arcs[n * i + j].as_slice();
            let l2 = logscores.arcs[n * j + i].as_slice();
            let s = log_sum_exp_two(l1, l2); // todo: Average for now
            //println!("Vals: i={} j={} l1len={} l2len={} s1={} s2={} s={}", i, j, l1.len(), l2.len(), s1, s2, s);
            merge.push((i as u32, j as u32, s));
        }
        let mut ranking = merge;
        ranking.par_sort_unstable_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        Ok(Ranking { entries: ranking })
    }

    /// Compute ranking from logscore matrix.
    ///
    /// todo: Using a shared variable iter over multiple threads seems to perform much worse than the single-threaded version
    ///
    /// *Implements steps: 3, 4, 5*
    pub fn ranking_from_logscore_matrix_par(logscores: LogscoreMatrix) -> Result<Ranking, String> {
        let n = logscores.n;
        let arc_iter = VariableIter::new(n, 2);
        let arc_iter = Mutex::new(arc_iter);
        let merge = vec![];
        let merge = Mutex::new(merge);
        let nthreads = thread::available_parallelism().unwrap().get();
        thread::scope(|s| {
            for _ in 1..nthreads {
                s.spawn(|| {
                    loop {
                        let mut lock = arc_iter.lock();
                        if !lock.next() {
                            break;
                        }
                        let i = lock.at()[0] as usize;
                        let j = lock.at()[1] as usize;
                        drop(lock);
                        let l1 = logscores.arcs[n * i + j].as_slice();
                        let l2 = logscores.arcs[n * j + i].as_slice();
                        let s1 = log_sum_exp(l1);
                        let s2 = log_sum_exp(l2);
                        let s = (s1 + s2) / 2.0; // todo
                        //println!("Vals: i={} j={} l1len={} l2len={} s1={} s2={} s={}", i, j, l1.len(), l2.len(), s1, s2, s);
                        let mut lock = merge.lock();
                        lock.push((i as u32, j as u32, s));
                        drop(lock);
                    }
                });
            }
        });
        let mut ranking = merge.into_inner();
        ranking.par_sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        Ok(Ranking { entries: ranking })
    }

    #[test]
    #[ignore]
    fn view_results() {
        let file = File::open("sample.100.mn2").unwrap();
        use std::time::Instant;
        let t0 = Instant::now();
        println!("Counting...");
        let (variables, counts) = count_map_from_file_par(file, 2).unwrap();
        let n = variables.len();
        let t1 = Instant::now();
        println!("Took {}", (t1 - t0).as_millis());
        println!("Accumulating logscores...");
        let (logscores, ..) = logscores_from_map(2, &variables, counts, Some(32.0)).unwrap();
        let t2 = Instant::now();
        println!("Took {}", (t2 - t1).as_millis());
        println!("Ranking...");
        let ranking = ranking_from_logscore_matrix(logscores).unwrap();
        let t3 = Instant::now();
        println!("Took {}", (t3 - t2).as_millis());
        println!("{} k={} @ {}ms count={}ms acc={}ms rank={}ms", "sample.100.mn2", 2, (t3 - t0).as_millis(), (t1 - t0).as_millis(), (t2 - t1).as_millis(), (t3 - t2).as_millis());
        let mut c = 0;
        for r in ranking.entries {
            println!("{} [{}->{}, {}]", c, r.0, r.1, r.2);
            c += 1;
            if c == 100 {
                break;
            }
        }
    }

    /// todo Iterate in a different order that can be taken advantage of by multiple threads.
    fn ranking_from_logscore_matrix_par2(_logscores: LogscoreMatrix) -> Result<Ranking, String> {
        todo!()
    }

    /// Benchmark discovery
    ///
    /// Run this test, with logging:
    /// cargo test bench_kn -- --ignored --nocapture
    #[test]
    #[ignore]
    fn bench_kn() {
        bench_kn_w("sample.2500.mn2", 1, true);
        bench_kn_w("sample.2500.mn2", 1, false);
        bench_kn_w("sample.2500.mn2", 2, true);
        bench_kn_w("sample.2500.mn2", 2, false);
    }

    #[cfg(test)]
    fn bench_kn_w(filen: &str, k: u32, par: bool) {
        use std::time::Instant;
        use std::hint::black_box;
        let file = File::open(filen).unwrap();
        let o = if par {
            let t0 = Instant::now();
            println!("Counting...");
            let (variables, counts) = count_map_from_file_par(file, k).unwrap();
            let n = variables.len();
            let t1 = Instant::now();
            println!("Took {}", (t1 - t0).as_millis());
            println!("Accumulating logscores...");
            let (logscores, ..) = logscores_from_map(k, &variables, counts, Some(32.0)).unwrap();
            let t2 = Instant::now();
            println!("Took {}", (t2 - t1).as_millis());
            println!("Ranking...");
            let ranking = ranking_from_logscore_matrix_par(logscores).unwrap();
            let t3 = Instant::now();
            println!("Took {}", (t3 - t2).as_millis());
            println!("{} k={} par={} @ {}ms count={}ms acc={}ms rank={}ms", filen, k, par, (t3 - t0).as_millis(), (t1 - t0).as_millis(), (t2 - t1).as_millis(), (t3 - t2).as_millis());
            ranking

        } else {
            let t0 = Instant::now();
            println!("Counting...");
            let (variables, counts) = count_map_from_file(file, k).unwrap();
            let n = variables.len();
            let t1 = Instant::now();
            println!("Took {}", (t1 - t0).as_millis());
            println!("Accumulating logscores...");
            let (logscores, ..) = logscores_from_map(k, &variables, counts, Some(32.0)).unwrap();
            let t2 = Instant::now();
            println!("Took {}", (t2 - t1).as_millis());
            println!("Ranking...");
            let ranking = ranking_from_logscore_matrix(logscores).unwrap();
            let t3 = Instant::now();
            println!("Took {}", (t3 - t2).as_millis());
            println!("{} k={} par={} @ {}ms count={}ms acc={}ms rank={}ms", filen, k, par, (t3 - t0).as_millis(), (t1 - t0).as_millis(), (t2 - t1).as_millis(), (t3 - t2).as_millis());
            ranking
        };
        black_box(o);
    }

}

/// Markov-chain Monte-Carlo based discovery methods
pub mod mcmc {

    pub fn rank_mcmc() {
        todo!()
    }

}
