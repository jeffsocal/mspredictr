use extendr_api::prelude::*;
use regex::Regex;

/// Return the mass of a peptide sequence
/// @export
#[extendr]
fn peptide_mass(s: &str) -> f32 {
    let re = Regex::new(r"([A-Z])(\-*\d*\.*\d*)").unwrap();
    let mut m: f32 = 0.0;
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

/// Return the length of a peptide sequence
/// @export
#[extendr]
fn peptide_len(s: &str) -> usize {
    let re = Regex::new(r"([A-Z])(\-*\d*\.*\d*)").unwrap();
    let cap: Vec<_> = re.captures_iter(&s).collect();
    return cap.len();
}

fn mass_letter(r: &char) -> f32 {
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
        _   =>   0.0,
    }
}

/// Return the mass of an atom
/// @export
#[extendr]
fn mass_atomic(r: &str) -> f32 {
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
fn mass_proton() -> f32 {
    return 1.00727646688;
}

/// Return the mass of a amino acid residue
/// @export
#[extendr]
fn mass_residue(s: &str) -> f32 {
    let mut m = 0.0;
    for ch in s.chars() {
        m += mass_letter(&ch);
    }
    return m;
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rmstandem;
    fn peptide_mass;
    fn mass_atomic;
    fn mass_proton;
    fn mass_residue;
    fn peptide_len;
}
