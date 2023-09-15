//! Iterators over variables, neighbours and configurations.

use std::ops::{Index, IndexMut};
use crate::{Array, Config, Variables};

/// An iterator over all `k`-combinations of `n` variables.
pub struct VariableIter {
    n: usize,
    at: Array<u32>,
    done: bool,
}

impl VariableIter {

    pub fn new(n: usize, k: usize) -> Self {
        if n == 0 || k == 0 || k > n {
            return Self { n, at: Box::new([]), done: true }
        }
        let mut variables = Vec::with_capacity(k);
        for i in 0..k-1 {
            variables.push(i as u32);
        };
        variables.push(k.wrapping_sub(2) as u32);
        let variables = variables.into_boxed_slice();
        VariableIter { n, at: variables, done: false }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn k(&self) -> usize {
        self.at.len()
    }

    pub fn at(&self) -> &[u32] {
        &self.at
    }

    /// Returns false when all combinations have been iterated through.
    pub fn next(&mut self) -> bool {
        if self.done {
            return false;
        }
        let mut pointer = self.k() - 1;
        loop {
            self.at[pointer] = self.at[pointer].wrapping_add(1);
            if self.at[pointer] as usize == self.n - self.k() + 1 + pointer {
                if pointer != 0 { // All digits are k.
                    pointer -= 1;
                    continue;
                } else {
                    self.done = true;
                    return false;
                };
            } else {
                if pointer == self.k() - 1 {
                    return true;
                } else {
                    pointer += 1;
                    self.at[pointer] = self.at[pointer - 1];
                    continue;
                };
            };
        };
    }

}

#[test]
fn test_var_iter() {
    test_iter(5, 3, 10, "[0 1 2][0 1 3][0 1 4][0 2 3][0 2 4][0 3 4][1 2 3][1 2 4][1 3 4][2 3 4]");
    test_iter(4, 3, 4, "[0 1 2][0 1 3][0 2 3][1 2 3]");
    test_iter(7, 1, 7, "[0][1][2][3][4][5][6]");
    test_iter(1, 1, 1, "[0]");
    test_iter(7, 7, 1, "[0 1 2 3 4 5 6]");
    test_iter(0, 3, 0, "");
    test_iter(3, 5, 0, "");
}

#[cfg(test)]
fn test_iter(n: usize, k: usize, vs: usize, eq: &str) {
    let mut iter = VariableIter::new(n, k);
    let mut c = 0;
    let mut s = String::new();
    while iter.next() {
        c += 1;
        let mut st = format!("[");
        for e in iter.at.iter() {
            st = format!("{}{} ", st, e);
        }
        st = format!("{}", st.trim());
        st = format!("{}]", st);
        s = format!("{}{}", s, st);
    }
    assert_eq!(c, vs);
    assert_eq!(s, eq)
}

/// Child config iterator
///
/// Given a parent config, iterate over all children configs.
pub struct ChildConfigIter<'a> {
    n: u32,
    variables: &'a Variables,
    buffer: Array<Config>,
    index: usize,
    child_cardinality: u32,
    done: bool,
}

impl IndexMut<usize> for ChildConfigIter<'_> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.buffer[index]
    }

}

impl Index<usize> for ChildConfigIter<'_> {

    type Output = Config;

    fn index(&self, index: usize) -> &Self::Output {
        &self.buffer[index]
    }

}

impl <'a> ChildConfigIter<'a> {

    pub fn new(variables: &'a Variables, parents: &'a [Config]) -> Self {
        let mut buffer = Vec::new();
        buffer.push(Config { variable: u32::MAX, value: 0 });
        for parent in parents {
            buffer.push(parent.clone());
        }
        let buffer = buffer.into_boxed_slice();
        Self { n: variables.len() as u32, variables, buffer, index: 0, child_cardinality: 0, done: false }
    }

    pub fn next(&mut self) -> bool {
        if self.done {
            return false;
        }
        let mut index = self.index;
        if self[index].value != self.child_cardinality { // There are more child configs to iterate over.
            self[index].value += 1;
            return true;
        }
        self[index].variable = self[index].variable.wrapping_add(1);
        loop {
            if index == self.buffer.len() - 1 { // Passed all parents.
                if self[index].variable == self.n { // All child variables have been iterated over.
                    self.done = true;
                    return false;
                } else { // This child variable has not been iterated over.
                    break;
                }
            } else {
                if self[index].variable == self[index + 1].variable { // Child is equal to a parent.
                    self[index] = self[index + 1];
                    index += 1;
                    self.index += 1;
                    self[index].variable += 1;
                } else {
                    break;
                }
            }
        }
        self[index].value = 0;
        self.child_cardinality = self.variables[self[index].variable];
        true
    }

    pub fn child(&self) -> Config {
        self[self.index]
    }

    pub fn child_cardinality(&self) -> u32 {
        self.child_cardinality
    }

    pub fn buffer(&self) -> &[Config] {
        &self.buffer
    }

    fn index(&self) -> usize {
        self.index
    }

}

//
// /// An iterator over all neighbour sets of x with cardinality k.
// ///
// /// # Fields
// /// * k: Number of neighbours.
// /// * n: Number of variables.
// /// * x: The neighbourhood variable.
// ///
// /// The initial configuration is [0, 1, ..., x - 1, x + 1, ..., m], with pointer beyond the last entry.
// ///
// pub struct NeighbourIterator {
//     k: u32,
//     n: u32,
//     x: u32,
//     configuration: Box<[u32]>,
// }
//
// impl NeighbourIterator {
//
//     /// Create a new NeighbourIterator.
//     ///
//     /// The
//     ///
//     /// Some combinations of `k`, `x`, `n` do not produce any valid configurations. In these
//     /// cases, [None] is returned.
//     pub fn new(k: u32, x: u32, n: u32) -> Option<Self> {
//         let mut configurations = Vec::with_capacity(k as usize);
//         let mut next = 0;
//         for _ in 0..k {
//             // Skip x.
//             if next == x {
//                 next += 1;
//             };
//             // Not enough variables to have k neighbours.
//             if next == n {
//                 return None;
//             };
//             configurations.push(next);
//             next += 1;
//         };
//         Some(Self { k, n, x, configuration: configurations.into_boxed_slice() })
//     }
//
//     /// Move to the next valid configuration. Returns false when all valid
//     /// configurations have been iterated over.
//     ///
//     /// Assumes that the current configuration is valid.
//     pub fn next(&mut self) -> bool {
//         let mut pointer = self.k - 1 as usize;
//         loop {
//             self.configuration[pointer] += 1;
//             if self.configuration[pointer] == self.n {
//                 // All digits are k.
//                 if pointer == 0 {
//                     return false;
//                 };
//                 pointer -= 1;
//                 continue;
//             } else if self.configuration[pointer] == self.x {
//                 // Always skip over x.
//                 continue;
//             } else {
//                 if pointer == self.k - 1 {
//                     return true;
//                 } else {
//                     pointer += 1;
//                     self.configuration[pointer] = self.configuration[pointer - 1];
//                     continue;
//                 };
//             };
//         };
//     }
//
//     pub fn print_config(&self) -> String {
//         let mut str = String::new();
//         str = format!("[");
//         for i in self.configuration {
//             str = format!("{}{} ", str, i);
//         };
//         str = format!("{}{} ", str, "]");
//         str
//     }
//
// }
//
// pub struct ConfigurationIter<'a> {
//     categories: &'a[u32],
//     configuration: Box<[u32]>,
// }
//
// impl <'a> ConfigurationIter<'a> {
//
//     pub fn next(&mut self) -> bool {
//         let n = self.configuration.len();
//         let mut pointer = n - 1;
//         loop {
//             self.configuration[pointer] += 1;
//             if self.configuration[pointer] == self.categories[pointer] {
//                 if pointer == 0 {
//                     return false;
//                 };
//                 pointer -= 1;
//                 continue;
//             };
//             if pointer == n - 1 {
//                 return true;
//             } else {
//                 pointer += 1;
//                 self.configuration[pointer] = u32::MAX; // Will wrap to 0 in next iteration.
//             };
//         };
//     }
//
//     /// Assumes that [0, 0, ..., 0] is a valid configuration.
//     pub fn new(variables: &[u32]) -> Self {
//         let configuration = vec![0; variables.len()].into_boxed_slice();
//         Self { categories: variables, configuration }
//     }
//
// }
