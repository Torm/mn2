//! Program for discovery of Markov networks. Behaviour is tuned by command line arguments.

use std::fs::File;
use std::str::FromStr;
use std::time::{Instant};
use console::Term;
use indicatif::{MultiProgress, ProgressBar};
use mn2::count::count_file_configs_into_map_par;
use mn2::discovery::k_neighbours::{logscores_from_map, ranking_from_logscore_matrix};

/// Use command: (can replace parameters 3 and sample.100.mn2)
/// cargo run --release -- 2 sample.100.mn2
fn main() {
    let mut args = std::env::args();
    args.next(); // The first arg is the binary. Skip.
    let k = match args.next() {
        None => {
            println!("No k value specified.");
            return;
        }
        Some(k) => u32::from_str(&k).unwrap(),
    };
    let path = match args.next() {
        None => {
            println!("No file path specified.");
            return;
        }
        Some(a) => a,
    };
    let file = File::open(&path).unwrap();
    let terminal = Term::stdout();
    terminal.write_line("#########################").unwrap();
    terminal.write_line(&format!("Discovery method=kn k={} file={}", k, &path)).unwrap();
    terminal.write_line("#########################").unwrap();
    terminal.write_line("[ ] Count...").unwrap();
    let bars = MultiProgress::new();
    let row_bar = bars.add(ProgressBar::new(1000));
    let t0 = Instant::now();
    let (variables, counts) = count_file_configs_into_map_par(file, k as usize + 1, Some((&bars, &row_bar))).unwrap();
    let n = variables.len();
    let t1 = Instant::now();
    bars.clear().unwrap();
    terminal.move_cursor_left(10000).unwrap();
    terminal.move_cursor_up(1).unwrap();
    terminal.write_line("[X] Count").unwrap();
    terminal.write_line(&format!("        [time={}ms] [unique-configs={}]", (t1 - t0).as_millis(), counts.len())).unwrap();
    terminal.write_line("[ ] Accumulate logscores...").unwrap();
    let t0 = Instant::now();
    let (logscores, max, min, accumulated, discarded) = logscores_from_map(k, &variables, counts, Some(10000.0)).unwrap();
    let t1 = Instant::now();
    terminal.move_cursor_left(10000).unwrap();
    terminal.move_cursor_up(1).unwrap();
    terminal.write_line("[X] Accumulate logscores").unwrap();
    terminal.write_line(&format!("        [time={}ms] [accumulated={}] [discarded={}] [total={}] [max={}] [min={}]", (t1 - t0).as_millis(), accumulated, discarded, accumulated + discarded, max, min)).unwrap();
    terminal.write_line("[ ] Compute scores, merge scores, rank arcs...").unwrap();
    let t0 = Instant::now();
    let ranking = ranking_from_logscore_matrix(logscores).unwrap();
    let t1 = Instant::now();
    terminal.move_cursor_left(10000).unwrap();
    terminal.move_cursor_up(1).unwrap();
    terminal.write_line("[X] Compute scores, merge scores, rank arcs").unwrap();
    terminal.write_line(&format!("        [time={}ms]", (t1 - t0).as_millis())).unwrap();
    terminal.write_line(&format!("Top 100 rankings:")).unwrap();
    let mut c = 0;
    for r in ranking.iter() {
        terminal.write_line(&format!("{} [{} -> {}] [{}]", c, r.0, r.1, r.2)).unwrap();
        c += 1;
        if c == 100 {
            break;
        }
    }
}
