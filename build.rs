use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const LOOKUP_PATH: &str = "./database/human_v_lookup.csv";

fn main() {
    let outpath = Path::new(&env::var("OUT_DIR").unwrap()).join("lookup.rs");
    let mut output = BufWriter::new(File::create(outpath).unwrap());
    let input = BufReader::new(File::open(LOOKUP_PATH).unwrap());

    let mut lookup_str: String = "{match (gene1, gene2) {".to_owned();

    for line_res in input.lines() {
        let line = line_res.unwrap();
        let mut split = line.split(",");
        let g1: &str = split.next().unwrap();
        let g2: &str = split.next().unwrap();
        let dist: &str = split.next().unwrap();
        lookup_str.push_str(&format!("(b\"{g1}\", b\"{g2}\") => {dist},\n"));
    }

    lookup_str.push_str("_ => panic!(\"Invalid gene(s)/allele(s) or comparison: {}, {}.\", str::from_utf8(gene1).unwrap(), str::from_utf8(gene2).unwrap())");
    lookup_str.push_str("}}");

    writeln!(
        &mut output,
        "use std::str;\npub fn total_distance(gene1: &[u8], gene2: &[u8]) -> u16 {}",
        lookup_str
    )
    .unwrap();
}
