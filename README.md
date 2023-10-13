# mn2

Discovery of Markov network structures.

## Cargo usage

Run with command: `cargo run -- <k> <mn2-file>`. Replace `<k>` with the number of neighbours
and `<mn2-file>` with the path to the file containing the data.

## Binary format

The input format describes tabular data with a header. The header contains information about
the number of columns and the ranges of the variables.

- The first 32 bits encode the number of columns. This number is known as `n`.
- The next `n` 32 bits encode the maximum value of the `n` variables.
- The remaining bits encode the rows of the data. The number of bits must be a multiple
  of `n` * 32.

**Example:**

This is a visualisation of the format. Each number represents 32-bits. Ignore all
whitespace.

```
5         # Header n
1 1 2 1 1 # Header variable ranges
1 0 0 1 0 # Data row 1
1 0 2 1 0 # Data row 2
0 1 0 0 0 # Data row 3
0 1 2 1 0 # Data row 4
...
```

## Samples

- `sample.10000.mn2`: 100*100 grid network, 10000 binary variables, 1000 rows, 19999 unique k1
- `sample.1000.mn2`: 20*50 grid network, 1000 binary variables, 1000 rows, 5000 unique k1
- `sample.100.mn2`: 10*10 grid network, 100 binary variables, 1000 rows, 200 unique k1, 20000 unique
  k2, 1311247 unique k3
- `sample.10.mn2`: 2*5 grid network, 10 binary variables, 1000 rows
- `sample.4.mn2`: 2*2 grid network, 4 binary variables, 1000 rows

## Benchmarks

- The count benchmarks shows that computation is CPU-bound and that multithreading
  speeds this up. Perhaps GPU-computing can be used here?
