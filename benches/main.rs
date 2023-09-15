use std::fs::File;
use std::thread::available_parallelism;
use std::time::Instant;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mn2::count;

criterion_group!(benches, benchmark);
criterion_main!(benches);

fn bench_(){}

fn benchmark(c: &mut Criterion) {
    use parking_lot::Mutex;
    use std::sync::Arc;
    println!("Benchmarking w. av. par: {}", available_parallelism().unwrap());
    {
        let t0 = Instant::now();
        let file = File::open("sample.10000.mn2").unwrap();
        let c = count::count_file_configs_into_map(file, 1).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=0 d=1000*10000 t=single @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.10000.mn2").unwrap();
        let c = count::count_file_configs_into_map_par(file, 1).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=0 d=1000*10000 t=multi @ {}", (t1 - t0).as_millis());
    }
    // These are too intensive.
    // {
    //     let t0 = Instant::now();
    //     let file = File::open("../sample.10000.mn2").unwrap();
    //     count::count_file_configs_into_map(file, 2).unwrap();
    //     let t1 = Instant::now();
    //     println!("count k=2 d=1000*10000 t=single @ {}", (t1 - t0).as_millis());
    // }
    // {
    //     let t0 = Instant::now();
    //     let file = File::open("../sample.10000.mn2").unwrap();
    //     count::count_file_configs_into_map_par(file, 2).unwrap();
    //     let t1 = Instant::now();
    //     println!("count k=2 d=1000*10000 t=multi @ {}", (t1 - t0).as_millis());
    // }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map(file, 1).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=0 d=1000*100 t=single @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map_par(file, 1).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=0 d=1000*100 t=multi @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map(file, 2).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=1 d=1000*100 t=single @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map_par(file, 2).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=1 d=1000*100 t=multi @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map(file, 3).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=2 d=1000*100 t=single @ {}", (t1 - t0).as_millis());
    }
    {
        let t0 = Instant::now();
        let file = File::open("sample.100.mn2").unwrap();
        let c = count::count_file_configs_into_map_par(file, 3).unwrap();
        black_box(c);
        let t1 = Instant::now();
        println!("count k=2 d=1000*100 t=multi @ {}", (t1 - t0).as_millis());
    }
    // Requires >4GB memory
    // {
    //     let t0 = Instant::now();
    //     let file = File::open("sample.100.mn2").unwrap();
    //     let c = count::count_file_configs_into_map_par(file, 4).unwrap();
    //     black_box(c);
    //     let t1 = Instant::now();
    //     println!("count k=4 d=1000*100 t=multi @ {}", (t1 - t0).as_millis());
    // }
}