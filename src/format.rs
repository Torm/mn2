//! Read and write mn2 binary format.

use crate::Variables;

fn write(variables: Variables, rows: Vec<u32>) -> Vec<u32> {
    let mut out = vec![];
    out.push(variables.n() as u32);
    for v in variables.variables().iter() {
        out.push(*v);
    }
    for
    out
}
