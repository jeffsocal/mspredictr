use extendr_api::prelude::*;
use regex::Regex;

fn peptide_mass_single(s: &str) -> f64 {
    let re = Regex::new(r"([A-Z])(\-*\d*\.*\d*)").unwrap();
    let mut m: f64 = 0.0;
    let cap: Vec<_> = re.captures_iter(&s).collect();
    for c in &cap {
        let ai = &c[1].chars().as_str();
        let am = &c[2].to_string();
        let mr = mass_residue(&ai);
        m += mr;
        if am != "" {
            let mm = &am.parse().unwrap();
            m += mm;
        }
    }

    return m + mass_atomic("O") + mass_atomic("H") * 2.0;
}

/// Return the mass of a peptide sequence
/// @export
#[extendr]
fn peptide_xleucine(s: &str) -> String {
    let out  = str::replace(s, "I", "L");
    return out;
}

/// Return the mass of a peptide sequence
/// @export
#[extendr]
fn peptide_mass(sequences: Strings) -> Vec<f64> {
   let masses = sequences.iter().map(|x| peptide_mass_single(&x)).collect();
   return masses;
}

fn peptide_length_single(s: &str) -> usize {
    let re = Regex::new(r"([A-Z])(\-*\d*\.*\d*)").unwrap();
    let cap: Vec<_> = re.captures_iter(&s).collect();
    let n = cap.iter().count();
    return n;
}

/// Return the mass of a peptide sequence
/// @export
#[extendr]
fn peptide_length(sequences: Strings) -> Vec<usize> {
   let masses = sequences.iter().map(|x| peptide_length_single(&x)).collect();
   return masses;
}

fn mass_letter(r: &char) -> f64 {
    match r {
        'A' =>  71.037_12,
        'R' => 156.101_1,
        'N' => 114.042_93,
        'D' => 115.026_94,
        'C' => 103.009_186,
        'E' => 129.042_59,
        'Q' => 128.058_58,
        'G' =>  57.021_465,
        'H' => 137.058_91,
        'I' => 113.084_06,
        'L' => 113.084_06,
        'K' => 128.094_96,
        'M' => 131.040_48,
        'F' => 147.068_42,
        'P' =>  97.052_765,
        'S' =>  87.032_03,
        'T' => 101.047_676,
        'W' => 186.079_32,
        'Y' => 163.063_32,
        'V' =>  99.068_41,
        'U' => 150.953_63,
        'O' => 237.147_73,
        _   => 113.956,
    }
}

/// Return the mass of an atom
/// @export
#[extendr]
fn mass_atomic(r: &str) -> f64 {
    match r {
        "H"  =>   1.007825032,
        "He" =>   4.00260325,
        "Li" =>   7.01600405,
        "Be" =>   9.01218214,
        "B"  =>  11.00930555,
        "C"  =>  12.0,
        "N"  =>  14.00307401,
        "O"  =>  15.99491462,
        "F"  =>  18.99840321,
        "Ne" =>  19.99244018,
        "Na" =>  22.98976967,
        "Mg" =>  23.9850419,
        "Al" =>  26.98153844,
        "Si" =>  27.97692653,
        "P"  =>  30.97376151,
        "S"  =>  31.97207069,
        "Cl" =>  34.96885271,
        "Ar" =>  39.96238312,
        "K"  =>  38.96370693,
        "Ca" =>  39.96259123,
        "Sc" =>  44.95591021,
        "Ti" =>  47.94794711,
        "V"  =>  50.94396371,
        "Cr" =>  51.94051192,
        "Mn" =>  54.93804961,
        "Fe" =>  55.93494212,
        "Co" =>  58.93320022,
        "Ni" =>  57.93534792,
        "Cu" =>  62.92960112,
        "Ga" =>  68.9255813,
        "As" =>  74.92159642,
        "Br" =>  78.91833762,
        "Kr" =>  83.9115073,
        "Rb" =>  84.91178933,
        "Sr" =>  87.90561432,
        "Y"  =>  88.90584793,
        "Zr" =>  89.90470372,
        "Nb" =>  92.90637752,
        "Rh" => 102.9055043,
        "Ag" => 106.9050936,
        "In" => 114.9038785,
        "Sb" => 120.903818,
        "I"  => 126.9044684,
        "Cs" => 132.9054473,
        "Ba" => 137.9052413,
        "La" => 138.9063483,
        "Ce" => 139.9054343,
        "Pr" => 140.9076483,
        "Eu" => 152.9212263,
        "Tb" => 158.9253433,
        "Ho" => 164.9303193,
        "Tm" => 168.9342113,
        "Lu" => 174.9407679,
        "Ta" => 180.9479963,
        "Re" => 186.9557508,
        "Ir" => 192.9629243,
        "Au" => 196.9665523,
        "Tl" => 204.9744123,
        "Pb" => 207.9766363,
        "Bi" => 208.9803833,
        "Th" => 232.0380504,
        "Pa" => 231.0358789,
        "U"  => 238.0507826,
        _    =>   0.0,
    }
}

/// Return the mass of a proton
/// @export
#[extendr]
fn mass_proton() -> f64 {
    return 1.00727646688;
}

/// Return the mass of a amino acid residue
/// @export
#[extendr]
fn mass_residue(s: &str) -> f64 {
    let mut m = 0.0;
    for ch in s.chars() {
        m += mass_letter(&ch);
    }
    return m;
}

/// Return the charged mass
/// @export
#[extendr]
fn mass_charged(mass: f64, z: i8) -> f64 {
    if z == 0 {
      return mass;
    }
    let z = z as f64;
    let m = (mass / z) + mass_proton();
    return m;
}

/// Return the neutral mass
/// @export
#[extendr]
fn mass_neutral(mz: f64, z: i8) -> f64 {
    if z == 0 {
      return mz;
    }
    let z = z as f64;
    let m = (mz - mass_proton()) * z;
    return m;
}

/// Return the mass ladder vector
/// @export
#[extendr]
pub fn mass_ladder(seq: &str) -> Vec<f64> {

    let rpoly = Regex::new(r"[A-Z]|\[.+?\]").unwrap();
    let rptms = Regex::new(r"[A-Z]|\-*\+*\d+\.*\d*").unwrap();

    // convert the REGEX to a string vector
    let mut seg_mass: Vec<f64> = vec![];
    for segment in rpoly.captures_iter(&seq) {
        let segment_this = segment.get(0).unwrap().as_str().to_string();

        if segment_this.len() == 1 {
            seg_mass.push(mass_letter(&segment_this.chars().next().expect("no amino acid found")));
        } else {
            let mut ptm_mass = 0.0;
            for ptms in rptms.captures_iter(&segment_this) {
                let ptms_this = ptms.get(0).unwrap().as_str().to_string();
                if ptms_this.len() == 1 {
                    ptm_mass += mass_letter(&ptms_this.chars().next().expect("no amino acid found"));
                } else {
                    ptm_mass += &ptms_this.as_str().parse().unwrap();
                    seg_mass.push(ptm_mass);
                }
            }
        }
    }
  return seg_mass;
}

/// Return the fragment mz vector
/// @export
#[extendr]
pub fn mass_fragments(seq: &str) -> Vec<f64> {

    let mass_h: f64 = mass_atomic("H");
    // let mass_n: f32 = mass_atomic("N");
    let mass_o: f64 = mass_atomic("O");
    // let mass_c: f32 = mass_atomic("C");
    let mass_water: f64 = mass_h * 2.0 + mass_o;
    // let mass_amine: f32 = mass_h * 3.0 + mass_n;

    let mut mz_frags = vec![];

    // push the precursor onto the struct
    let pep_ladder = mass_ladder(&seq);
    let pep_n: Vec<_> = (0..(pep_ladder.len() - 1)).collect();

    for i in pep_n.iter() {

        // create the xyz ion series
        let ion_y: f64 = pep_ladder.iter().rev().take(*i+1).sum();

        // create the abc ion series
        let ion_b: f64 = pep_ladder.iter().take(*i+1).sum();

        mz_frags.push((ion_y + mass_water) / 1.0 + mass_proton());
        if i > &0 {
          mz_frags.push((ion_b) / 1.0 + mass_proton());
        }

    }

    return mz_frags;
}

/// Return the fragment mz vector
/// @export
#[extendr]
pub fn index_fragments(seq: &str, tolerance: f64) -> Vec<i64> {
    let i_frag = mass_fragments(&seq)
                 .iter()
                 .map(|x| indexr(x, tolerance))
                 .collect::<Vec<_>>();
    return i_frag;
}

/// Return a mass index.
/// @export
#[extendr]
fn index_mass(masses: Vec<f64>, tolerance: f64) -> Vec<i64> {
  let indexes:Vec<i64> = masses
                         .iter()
                         .map(|x| indexr(x, tolerance))
                         .collect();
  return indexes;
}

fn indexr(mass: &f64, tolerance: f64) -> i64 {
  // let m_int = 0.98;
  // let m_slp = 0.999_521_6;
  // let out = mass * m_slp + m_int;
  let out = mass / ( tolerance / 3.0 );
  return out.round() as i64;
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rmstandem;
    fn mass_atomic;
    fn mass_proton;
    fn mass_residue;
    fn mass_charged;
    fn mass_neutral;

    fn peptide_mass;
    fn peptide_length;

    fn mass_ladder;
    fn mass_fragments;

    fn index_fragments;
    fn index_mass;

    fn peptide_xleucine;
}
