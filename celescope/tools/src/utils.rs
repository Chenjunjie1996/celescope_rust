#![allow(unused, warnings)]
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Read};
use std::fs::File;
use std::error::Error;
use std::ops::Range;
use bam;
use bam::header::{Header, HeaderEntry};
use bam::record::RECORD_PAIRED;
use bio::alphabets;
use bio::io::bed::Reader;
use bio::io::fastq;
use bio::io::fasta::{self, Record};
use csv;
use std::process::Command;

const BARCODE_FILE_NAME: &str = "barcode.tsv";
const FILTERED_MATRIX_DIR_SUFFIX: (&str, &str) = ("filtered_feature_bc_matrix", "matrix_10X");

struct Samtools{
    in_bam: String,
    out_bam: String,
    threads: i8,
    temp_sam_file: String,
    debug: bool,
}

impl Samtools {
    fn new(in_bam: String, out_bam:String, threads:i8, temp_sam_file:String, debug: bool) -> Self {
        Self {
            in_bam,
            out_bam,
            threads,
            temp_sam_file,
            debug,
        }
    }

    fn samtools_sort(&self) {
        Command::new(cmd);

        let cmd = format!(
            "samtools sort {} -o {} --threads {}", self.in_bam, self.out_bam, self.threads
        );
        subprocess::Exec::shell(cmd);
    }

    fn samtools_index(&self) {
        let cmd = format!(
            "samtools index {}", self.in_bam
        );
        subprocess::Exec::shell(cmd);
    }

    fn sort_bam(&self) {
        let cmd = format!(
            "samtools index {}", self.out_bam
        );
        subprocess::Exec::shell(cmd);
    }

}

fn generic_open(file_name: String) -> Result<File, Box<dyn Error>> {
    let mut file_obj = File::open(&file_name)?;
    Ok(file_obj)
}

fn read_one_col(file: String) -> (Vec<String>, usize) {
    let mut df = std::fs::File::open(&file).unwrap();
    let mut contents = String::new();
    df.read_to_string(&mut contents).unwrap();

    let mut bc_vec = Vec::new();
    for i in contents.split('\n') {
        bc_vec.push(i.to_string());
    }
    bc_vec.pop();
    (bc_vec.clone(), bc_vec.len())
}

fn read_one_col_fl(file: String) -> Vec<String>{
    let mut df = std::fs::File::open(&file).unwrap();
    let mut contents = String::new();
    df.read_to_string(&mut contents).unwrap();

    let mut bc_vec = Vec::new();
    for i in contents.split(',') {
        if i.contains("\r"){
            let barcode = {
                if i.contains("exact_subclonotype_id"){
                    i.replace("exact_subclonotype_id\r\n", "")
                } else {
                    i.replace("1\r\n", "")
                }
            };
            bc_vec.push(barcode); 
        }
    }
    bc_vec.pop();
    bc_vec
}

fn read_fasta(fasta_file: String) -> std::collections::HashMap<String, String> {
    let mut hm = HashMap::new();

    let mut reader = fasta::Reader::from_file(&fasta_file).unwrap();
    // let mut result_vec = Vec::new();
    for result in reader.records(){
        let record = result.expect("Error during fasta record parsing");

        let head_msgs = record.id();
        let head_msgs_vec: Vec<&str> = record.id().split('_').collect();
        let seq = String::from_utf8_lossy(record.seq());
        let barcode = head_msgs_vec[0];
        let umi = head_msgs_vec[1];
        hm.insert(head_msgs.to_string(), seq.to_string());
    }
    hm
}

fn haming_distance(str1: String, str2: String) -> usize{
    let mut distance = 0;
    let (length, length2) = (str1.len(), str2.len());
    if length != length2 {
        panic!("string1 {} and string2 {} do not have same length", length, length2);
    }
    let str1_vec = str1.chars().collect::<Vec<_>>();
    let str2_vec = str2.chars().collect::<Vec<_>>();
    for i in 0..length{
        if str1_vec[i] != str2_vec[i] { distance+=1 }
    }
    distance
}

fn haming_correct(str1: String, str2: String) -> bool{
    let threshold = str1.len() / 10 + 1;
    if haming_distance(str1, str2) < threshold{
        return true;
    }
    false
}

fn format_number(number: usize) -> String {
    number.to_string()
}

fn get_matrix_file_path(matrix_dir:String, file_name:String) -> String {
    let file_path_vec = vec![
        format!("{}/{}", &matrix_dir, &file_name),
        format!("{}/{}.gz", &matrix_dir, &file_name),
        ];
    for i in &file_path_vec{
        if check_file(i) == 1 {
            return i.to_string()
        }
    }
    panic!("Not Found")
}

fn check_file(bible_path: &str) -> i32 {
    let f = File::open(bible_path); 
    let result = match f {
        Ok(file) => 1,
        Err(err) => 0 
    }; 
    result
}

fn get_barcode_from_matrix_dir(matrix_dir: String) -> (Vec<std::string::String>, usize) {
    let match_barcode_file = get_matrix_file_path(matrix_dir, (&BARCODE_FILE_NAME).to_string());
    let (match_barcode, n_match_barcode) = read_one_col(match_barcode_file);
    (match_barcode, n_match_barcode)
}

fn get_matrix_dir_from_match_dir(match_dir: String) -> Vec<String> {
    let mut matrix_dir_pattern_list = Vec::new();
    matrix_dir_pattern_list.push(format!("{}/05.count/{}", &match_dir, FILTERED_MATRIX_DIR_SUFFIX.0));
    matrix_dir_pattern_list.push(format!("{}/05.count/{}", &match_dir, FILTERED_MATRIX_DIR_SUFFIX.1));
    matrix_dir_pattern_list
}

fn read_json(file_name: &str) -> serde_json::Value{
    let f = File::open(file_name).unwrap();
    let values:serde_json::Value = serde_json::from_reader(f).unwrap();
    values
}

fn get_fasta_read_number(file_name: &str) -> u32 {
    let mut nb_reads = 0;
    let mut reader = fasta::Reader::from_file(file_name).unwrap();
    for result in reader.records(){
        let record = result.expect("Error during fasta record parsing");
        nb_reads += 1;
    }
    nb_reads
}

fn get_fastq_read_number(file_name: &str) -> u32 {
    let mut nb_reads = 0;
    let mut reader = fastq::Reader::from_file(file_name).unwrap();
    for result in reader.records(){
        let record = result.expect("Error during fasta record parsing");
        nb_reads += 1;
    }
    nb_reads
}

fn fastq_line(name:&str, seq:&str, qual:&str) -> String {
    format!("@{}\n{}\n+\n{}\n", name, seq, qual)
}

fn get_assay_text(assay: &String) -> String{
    "Single-cell".to_string() + assay
}

fn add_tag(in_bam: &str){
    let reader = bam::BamReader::from_path(in_bam, 4).unwrap();
    for record in reader {
        let record = record.unwrap();
        // let attr = HeaderEntry::entry_name(&record);
        // record.tags_mut().push_num(b"CB", value)
    }
}

fn samtools_sort(in_file: &str, out_file: &str) {
    let mut cmd = format!("samtools sort {} -o {} --thread 1", in_file, out_file);
    subprocess::Exec::shell(cmd);
}

fn reverse_complement(seq: &str) -> String{
    // result vector
    let mut rc_seq: String = String::with_capacity(seq.len());

    // iterate through the input &str
    for c in seq.chars().rev() {
        if c == 'A' {
            rc_seq.push('T');
        } else if c == 'C' {
            rc_seq.push('G');
        } else if c == 'G' {
            rc_seq.push('C');
        } else if c == 'T' {
            rc_seq.push('A');
        } else {
            rc_seq.push('X');
        }
    }
    rc_seq
}
fn main() {
    let text = b"ACAGCTCGATCGGTA$";
    println!("{:?}", text);

    let seq = "ATCG";
    let rev_seq = reverse_complement(seq);
    println!("{:?}",rev_seq);

    // C:/Users/admin/Desktop/bd_flvdj/cr/0916TAA-H_220919N_S91/filtered_contig.fasta
    let path = "C:/Users/admin/Desktop/bd_flvdj/cr/0916TAA-H_220919N_S91/filtered_contig.fasta";
    let mut reader = fasta::Reader::from_file(&path).unwrap();
    let mut nb_reads = 0;

    let rna_barcode = vec!["AGCCCAATCCACCTACCCACTACATAT", "TAGCTTGTAGCTAACTCAGTCTTGGCA", "AAGCTCAACGCTAACTCAGTCTTGGCA", "ACCCTGACCGCTAACTCAGTCTTGGCA"];

    // let mut result_vec = Vec::new();
    for result in reader.records(){
        let record = result.expect("Error during fasta record parsing");
        nb_reads += 1;

        let head_msgs: Vec<&str> = record.id().split('_').collect();
        let seq = String::from_utf8_lossy(record.seq());
        let barcode = head_msgs[0];
        let umi = head_msgs[1];

        if rna_barcode.contains(&barcode){
            println!("barcode: {}",barcode);
            println!("umi: {}", umi);           
        }
    }
    println!("{}", nb_reads);

    let json_path = "C:/Users/admin/Desktop/python/test.json";
    let values = read_json(json_path);
    let tag1 = &values["TAG_1"];
    let tag2 = &values["TAG_2"];
    println!("{:?}", tag1);
    println!("{:?}", tag2);




}
