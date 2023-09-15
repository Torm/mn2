//! Compute configuration counts found in data.

use std::fs::File;
use std::io::{BufReader, ErrorKind};
use std::sync::Mutex;
use std::thread::available_parallelism;
use std::time::Duration;
use byteorder::{BigEndian, ReadBytesExt};
use dashmap::DashMap;
use fxhash::{FxBuildHasher};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use crate::iter::VariableIter;
use crate::{Array, Config, Variables};

/// Read a file and count the configs.
pub fn count_file_configs_into_map(file: File, k: usize) -> Result<(Variables, DashMap<Array<Config>, u32, FxBuildHasher>), String> {
    let mut reader = BufReader::new(file);
    //let mut c = 0;
    let n = match reader.read_u32::<BigEndian>() {
        Ok(n) => n as usize,
        Err(e) => return Err(format!("Error reading header n: {}", e.to_string())),
    };
    let mut variables = vec![0; n];
    if let Err(e) = reader.read_u32_into::<BigEndian>(&mut variables) {
        return Err(format!("Error reading variables: {}", e.to_string()));
    }
    let variables = variables.into_boxed_slice();
    let counts = DashMap::<Array<Config>, u32, FxBuildHasher>::default();
    let mut key_buffer = vec![Config { variable: 0, value: 0 }; k].into_boxed_slice();
    let mut value_buffer = vec![0; n].into_boxed_slice();
    loop {
        // c += 1;
        // if c % 1 == 0 {
        //     println!("{} entries read!", c);
        // }
        if let Err(e) = reader.read_u32_into::<BigEndian>(&mut value_buffer) {
            match e.kind() {
                ErrorKind::UnexpectedEof => break,
                e => return Err(e.to_string()),
            }
        }
        for k in 1..k+1 {
            let mut variable_iter = VariableIter::new(n, k);
            while variable_iter.next() {
                for i in 0..k {
                    let variable = variable_iter.at()[i];
                    let value = value_buffer[variable as usize];
                    let config = Config { variable, value };
                    key_buffer[i] = config;
                }
                if let Some(mut d) = counts.get_mut(&key_buffer[0..k]) {
                    *d += 1;
                } else {
                    let mut key = vec![Config { variable: 0, value: 0 }; k];
                    key_buffer[0..k].clone_into(&mut key);
                    counts.insert(key.into_boxed_slice(), 1);
                }
            }
        }
    }
    println!("Unique configs: {}", counts.len());//////////////
    let variables = Variables { variables };
    Ok((variables, counts))
}

/// Same as [count_file_configs_into_map], but processes rows in parallel.
///
/// Possible libraries:
/// - unordered_dense
/// - chashmap
/// - flurry
/// - dashmap
pub fn count_file_configs_into_map_par(file: File, k: usize, bars: Option<(&MultiProgress, &ProgressBar)>) -> Result<(Variables, DashMap<Array<Config>, u32, FxBuildHasher>), String> {
    let mut reader = BufReader::new(file);
    let n = match reader.read_u32::<BigEndian>() {
        Ok(n) => n as usize,
        Err(e) => return Err(format!("Error reading n: {}", e.to_string())),
    };
    let mut variables = vec![0; n];
    if let Err(e) = reader.read_u32_into::<BigEndian>(&mut variables) {
        return Err(format!("Error reading variables: {}", e.to_string()));
    }
    let variables = variables.into_boxed_slice();
    let counts = DashMap::<Array<Config>, u32, FxBuildHasher>::default();
    let tc = available_parallelism().unwrap().into();
    let mut c = 0;
    let reader = Mutex::new((reader, &mut c));
    std::thread::scope(|s| {
        for _ in 0..tc {
            let bars = if let Some((config_bars, row_bar)) = bars {
                let cfgs = (n-k+1..n+1).product::<usize>() / (1..k+1).product::<usize>();
                let config_bar = config_bars.add(ProgressBar::new(cfgs as u64).with_style(ProgressStyle::with_template("{spinner} {percent}%").unwrap()));
                config_bar.enable_steady_tick(Duration::from_millis(250));
                Some((row_bar, config_bar))
            } else {
                None
            };
            s.spawn(|| read_stream_par_task(n, k, &reader, &counts, bars));
        }
    });
    let variables = Variables { variables };
    Ok((variables, counts))
}

fn read_stream_par_task(n: usize, k: usize, reader: &Mutex<(BufReader<File>, &mut usize)>, counts: &DashMap<Array<Config>, u32, FxBuildHasher>, bar: Option<(&ProgressBar, ProgressBar)>) -> Result<(), String> {
    let mut key_buffer = vec![Config { variable: 0, value: 0 }; k].into_boxed_slice();
    let mut value_buffer = vec![0; n].into_boxed_slice();
    loop {
        let mut reader = reader.lock().unwrap();
        if let Err(e) = reader.0.read_u32_into::<BigEndian>(&mut value_buffer) {
            match e.kind() {
                ErrorKind::UnexpectedEof => break,
                e => return Err(e.to_string()),
            }
        }
        let mut c: u64 = 0; // Number of configs read from this row.
        if let Some((_row_bar, config_bar)) = &bar {
            config_bar.reset();
        }
        drop(reader);
        for k in 1..k+1 {
            let mut variable_iter = VariableIter::new(n, k);
            while variable_iter.next() {
                for i in 0..k {
                    let variable = variable_iter.at()[i];
                    let value = value_buffer[variable as usize];
                    let config = Config { variable, value };
                    key_buffer[i] = config;
                }
                if let Some(mut d) = counts.get_mut(&key_buffer[0..k]) {
                    *d += 1;
                } else {
                    let mut key = vec![Config { variable: 0, value: 0 }; k];
                    key_buffer[0..k].clone_into(&mut key);
                    counts.insert(key.into_boxed_slice(), 1);
                }
                c += 1;
                if c % 10000 == 0 {
                    if let Some((_row_bar, config_bar)) = &bar {
                        config_bar.inc(10000);
                    }
                }
            }
        }
        if let Some((row_bar, _config_bar)) = &bar {
            row_bar.inc(1);
        }
    }
    Ok(())
}

pub fn count_file_configs_into_tree(_file: File, _k: u32) {
    todo!()
}

pub fn count_file_configs_into_tree_par(_file: File, _k: u32) {
    todo!()
}

// /// A structure for the accumulation of counts.
// pub trait CountAccumulation {
//     fn increment(&mut self, config: &[Config]);
// }
//
// pub struct HashedCounts(DashMap<Box<[Config]>, u64>);
//
// impl CountAccumulation for HashedCounts {
//     fn increment(&mut self, _config: &[Config]) {
//
//     }
// }
//
// /// A collection of counts that can be iterated over.
// pub trait Counts {
// }
