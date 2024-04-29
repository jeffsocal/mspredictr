use extendr_api::prelude::*;
use regex::Regex;

mod isotopes;
mod cosine_similarity;

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

/// Return a I/L substituted sequence
/// @param s
/// The peptide sequence string.
/// @export
/// @examples
/// peptide_xleucine('SILLY')
#[extendr]
fn peptide_xleucine(s: &str) -> String {
    let out  = str::replace(s, "I", "L");
    return out;
}

/// Return the mass of a peptide sequence
/// @param sequences
/// A vector of peptide sequence strings.
/// @export
/// @examples
/// peptide_mass(c('SAMPLE', 'SILLY'))
#[extendr]
fn peptide_mass(sequences: Strings) -> Vec<f64> {
   let masses = sequences.iter().map(|x| peptide_mass_single(&x)).collect();
   return masses;
}

/// Return the peptide sequence reversed
/// @param sequences
/// A vector of peptide sequence strings.
/// @export
/// @examples
/// peptide_mass(c('THINK', 'SAM[15.99]PLER'))
#[extendr]
fn peptide_reverse(sequences: Strings) -> Vec<String> {
   let out = sequences.iter().map(|x| peptide_reverse_single(&x)).collect();
   return out;
}

fn peptide_reverse_single(seq: &str) -> String {

    let peptide = Regex::new(r"[A-Z]|\[.+?\]").unwrap();
    // convert the REGEX to a string vector
    let mut residues: Vec<_> = peptide
                            .captures_iter(&seq)
                            .map(|a| a[0].to_string())
                            .collect();

    // Reverse the vector
    residues.reverse();
    let residue_first = residues.get(0).unwrap();

    residues.push(residue_first.to_string());
    residues.remove(0);
    let str_peptide = residues.into_iter().collect::<String>();

  return str_peptide;
}

/// Return the unit length of a peptide sequence
/// @param sequences
/// A vector of peptide sequence strings.
/// @export
/// @examples
/// peptide_mass(c('SA[M15.99]PLE', 'SAM[15.99]PLE'))
#[extendr]
fn peptide_length(sequences: Strings) -> Vec<usize> {
   let masses = sequences.iter().map(|x| peptide_length_single(&x)).collect();
   return masses;
}

fn peptide_length_single(s: &str) -> usize {
    let re = Regex::new(r"[A-Z]|\[.+?\]").unwrap();
    let cap: Vec<_> = re.captures_iter(&s).collect();
    let n = cap.iter().count();
    return n;
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
/// @param atom
/// The IUPAC amomic letter designation
/// @export
/// @examples
/// mass_atomic('C')
/// mass_atomic('Na')
#[extendr]
fn mass_atomic(atom: &str) -> f64 {
    match atom {
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
/// @examples
/// mass_proton()
#[extendr]
fn mass_proton() -> f64 {
    return 1.00727646688;
}

/// Return the mass of a neutron
/// @export
/// @examples
/// mass_neutron()
#[extendr]
fn mass_neutron() -> f64 {
  return 1.00866491588;
}

/// Return the mass of a amino acid residue
/// @param s
/// The single letter designation for an amino acid
/// @export
/// @examples
/// mass_residue('G')
/// mass_residue('Y')
#[extendr]
fn mass_residue(s: &str) -> f64 {
    let mut m = 0.0;
    for ch in s.chars() {
        m += mass_letter(&ch);
    }
    return m;
}

/// Return the charged mass
/// @param mass
/// The floating point mass value
/// @param z
/// The integer charge value
/// @export
/// @examples
/// mass_charged(1234.567, 2)
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
/// @param mz
/// The floating point mass value
/// @param z
/// The integer charge value
/// @export
/// @examples
/// mass_neutral(1234.567, 2)
#[extendr]
fn mass_neutral(mz: f64, z: i8) -> f64 {
    if z == 0 {
      return mz;
    }
    let z = z as f64;
    let m = (mz - mass_proton()) * z;
    return m;
}

/// Return the vector of mass values for each peptide unit
/// @param seq
/// The peptide sequence string
/// @export
/// @examples
/// mass_ladder('SA[M15.99]PLER')
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

/// Return a simplified fragment mz vector
/// @param seq
/// The peptide sequence string
/// @export
/// @examples
/// mass_fragments('SA[M15.99]PLER')
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
        mz_frags.push((ion_b) / 1.0 + mass_proton());

    }

    return mz_frags;
}

/// Return the fragment indexes
/// @param seq
/// The peptide sequence string
/// @param tolerance
/// The numerical float for mass tolerance in Th
/// @export
/// @examples
/// index_fragments('SA[M15.99]PLER', 0.05)
#[extendr]
pub fn index_fragments(seq: &str, tolerance: f64) -> Vec<i64> {
    let i_frag = mass_fragments(&seq)
                 .iter()
                 .map(|x| (x / tolerance) as i64)
                 .collect::<Vec<_>>();
    return i_frag;
}

/// Return a mass indexes.
/// @param masses
/// A vector of mass sequences
/// @param tolerance
/// The numerical float for mass tolerance in Th
/// @export
/// @examples
/// index_mass(c(123.45, 234.56, 345.67), 0.05)
#[extendr]
fn index_mass(masses: Vec<f64>, tolerance: f64) -> Vec<i64> {
  let indexes:Vec<i64> = masses
                         .iter()
                         .map(|x| (x / tolerance) as i64)
                         .collect();
  return indexes;
}

/// Returns the boolean index of the top N largest values
/// @param f
/// A vector of numerical floats
/// @param n
/// The number of top values to keep
/// @export
/// @examples
/// which_top_n(c(123.45, 234.56, 345.67), 2)
#[extendr]
fn which_top_n(f: Vec<f64>, n: i32) -> Vec<bool> {
    let fl = f.len();
    if fl <= n as usize {
      let ft = f.iter().map(|_x| 1 == 1 ).collect::<Vec<_>>();
      return ft;
    }

    let mut fs = f.clone();
    // sort
    fs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // find the n-th value
    let fi = fs[fl - n as usize];
    let ft = f.iter().map(|x| x >= &fi ).collect::<Vec<_>>();

    return ft.to_vec();
}

/// Returns the boolean index of values that are in proximity to mz
/// @param f
/// A vector of numerical floats
/// @param mz
/// The the value to look for
/// @export
/// @examples
/// which_xprecursor(c(123.45, 234.56, 345.67), 233.61)
#[extendr]
fn which_xprecursor(f: Vec<f64>, mz: f64) -> Vec<bool> {
    let g = f.iter().map(|x| (x - mz).abs() > 2.5 ).collect::<Vec<_>>();
    return g;
}

use crate::isotopes::isotopes;
use crate::cosine_similarity::cosine_similarity;

// Take mz features and derrive isotope clusters to form molecular fetaures
//
// 1. take the tallest unassigned peak
// 2. find all peaks within the mz & rt window
// 3. score against the isotope profiles from different charge states
// 4. assign to the best fit
// 5. repeat 1 until all peaks are assigned

#[extendr]
struct IsotopeFeature {
    pub z: i8,
    pub mz: f64,
    pub intensity: f64,
    pub cluster: i64,
    pub isoscore: f64
}

fn build_point(mz: f64, intensity: f64) -> IsotopeFeature {
    IsotopeFeature {
        z: 0,
        mz,
        intensity,
        cluster: 0,
        isoscore: 0.0
    }
}

#[derive(Debug)]
struct IsotopicFit {
    indexes: Vec<usize>,
    charge: i8,
    score: f64,
}


/// Returns the isotopic grouping for a given mass spectrum take in only mz, rt,
/// and abundance as vectors then build the Vec<IsotopeFeature> this gets around
/// the use of data.frame in R
/// @param vec_mz
/// A vector of numerical floats representing the mz component (both vectors must be sorted on this value)
/// @param vec_int
/// A vector of numerical floats representing the ion intensity component
/// @export
/// @examples
/// group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137, 291.164, 292.177, 293.124, 296.135, 298.139),
///                c(218487, 44736, 29195, 1021168, 46029, 104552, 21997, 15262, 19908, 61741))
#[extendr]
fn group_isotopes(vec_mz: Vec<f64>,
                  vec_int: Vec<f64>) -> Vec<i64> {

    let mut points: Vec<IsotopeFeature> = vec_mz
        .into_iter()
        .zip(vec_int)
        .map(|(mz, intensity)| build_point (
            mz,
            intensity
        ))
        .collect();

    // iterate through all the points
    for i in 0..points.len() {

        if points[i].cluster != 0 {
            continue;
        };

        let this_mz = points[i].mz;
        let this_intensity = points[i].intensity;

        let mut isotope_fits: Vec<IsotopicFit> = vec![];

        // iterate over each charge state
        for charge in vec![6,5,4,3,2,1]{

            // get the expected isotopic profile
            let isotopes = isotopes(this_mz, charge, this_intensity);

            let mut matches: Vec<usize> = vec![];
            let mut intensity_observe: Vec<f64> = vec![];
            let mut intensity_predict: Vec<f64> = vec![];

            // iterate over each isotope
            for isotope in &isotopes {

                // keep track of the matching isotopes
                // if an isotope is not matched drop out of the loop
                let mut has_match = false;

                // iterate over every point in the vector
                for ii in i..points.len() {
                    let x = &points[ii];

                    if x.mz < isotope.mz - 0.1 {continue;}
                    if x.mz > isotope.mz + 0.1 {break;}

                    if x.cluster == 0
                        && x.intensity < isotope.intensity * 3.33
                        && x.intensity > isotope.intensity * 0.33 {

                        // collect the indexes that match isotopes
                        matches.push(ii);

                        // collect the abundances for regression analysis
                        intensity_observe.push(points[ii].intensity);
                        intensity_predict.push(isotope.intensity);

                        // register the match
                        has_match = true;

                        // stop iterating on points, go to next isotope
                        break;
                    }
                }

                // if an isotope is not matched drop out of the loop
                if has_match == false {
                    break;
                }
            }

            // perform a regression on any clusters
            let match_len: f64 = matches.len() as f64;
            if match_len == 1.0 {
                isotope_fits.push(
                    IsotopicFit{indexes: matches,
                        charge,
                        score: 0.0}
                )
            } else if match_len > 1.0 {
                // perform a linear regression on the abundance values
                // to determine best fit for the isotopic distribution
                // let fit = linear_regression(abundnce_predict, abundnce_observe);
                let fit = cosine_similarity(intensity_predict, intensity_observe);
                isotope_fits.push(
                    IsotopicFit{indexes: matches,
                        charge,
                        score: fit * match_len}
                )
            }
        }

        // if the feature has any matching features take the
        // best fitting isotope pattern
        if isotope_fits.len() > 0 {
            isotope_fits.sort_by_key(|i| ordered_float::OrderedFloat(-i.score));
            for m in &isotope_fits[0].indexes {
                points[*m].isoscore = (isotope_fits[0].score).round();
                points[*m].cluster = i as i64 + 1;
                points[*m].z = isotope_fits[0].charge as i8;
            }
        }
    }

    //points.sort_by_key(|p| ordered_float::OrderedFloat(p.cluster as f64 * p.mz));

    points.iter().map(|x| x.cluster).collect()
}

/// A Helper function that returns the array index of the monoisotopes given group_isotopes()
/// @param vec_iso
/// A vector of numerical floats representing the isotopic groups
/// @export
/// @examples
/// group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137, 291.164, 292.177, 293.124, 296.135, 298.139),
///                c(218487, 44736, 29195, 1021168, 46029, 104552, 21997, 15262, 19908, 61741)) |>
///   which_monoisotopes()
#[extendr]
fn which_monoisotopes(vec_iso: Vec<f64>) -> Vec<i64> {

    let mut points: Vec<i64> = vec![];
    let mut isotopes: Vec<f64> = vec![];
    for i in 0..vec_iso.len() {
      if isotopes.contains(&vec_iso[i] as &f64)  {
        continue;
      }
      isotopes.push(vec_iso[i]);
      points.push(i as i64 + 1);
    }
    return points;
}


/// Helper function that provides the isotopic assignment given the output from group_isotopes()
/// @param vec_iso
/// A vector of numerical floats representing the isotopic groups
/// @export
/// @examples
/// group_isotopes(c(287.171, 288.119, 288.174, 290.161, 291.137, 291.164, 292.177, 293.124, 296.135, 298.139),
///                c(218487, 44736, 29195, 1021168, 46029, 104552, 21997, 15262, 19908, 61741)) |>
///   label_isotopes()
#[extendr]
fn label_isotopes(vec_iso: Vec<f64>) -> Vec<f64> {

    let mut points: Vec<f64> = vec![];
    for i in 0..vec_iso.len() {
      let iso_id = vec_iso[i] as f64;
      let iso_no = vec_iso[0..i].iter().filter(|&n| *n == iso_id).count();
      points.push(iso_no as f64);
    }
    return points;
}



// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod mspredictr;
    fn mass_atomic;
    fn mass_proton;
    fn mass_neutron;
    fn mass_residue;
    fn mass_charged;
    fn mass_neutral;

    fn peptide_mass;
    fn peptide_length;
    fn peptide_reverse;

    fn mass_ladder;
    fn mass_fragments;

    fn index_fragments;
    fn index_mass;

    fn peptide_xleucine;

    fn which_top_n;
    fn which_xprecursor;
    fn which_monoisotopes;

    fn group_isotopes;
    fn label_isotopes;
}
