#' helper function to read in platform specific results
#'
#' @param x location of file to parse
#'
#' @return a tibble
#'
read_omssa <- function(
    x,
    cpus = 1
){

  proton_mass <- mass_proton()

  # Filename/id contains experiment scan level (eg ms1, ms2 included)
  out <- x |> readr::read_csv(show_col_types = FALSE) |>
    dplyr::mutate(psm_dp = NA,
                  psm_peptide = NA,
                  psm_score = `E-value` |> log10() * -1) |>
    dplyr::rename(
      # 1Th correction to get [M+H]+
      psm_mh = Mass + proton_mass,
      ms_event = `Filename/id`,
      psm_sequence = Peptide,
      psm_protein = Defline
    ) |>
    dplyr::mutate(ms_event = stringr::str_remove(ms_event, ".+scan\\=") |> as.numeric(),
                  psm_protein = psm_protein |> stringr::str_remove("\\s.+")) |>
    dplyr::group_by(ms_event) |>
    dplyr::mutate(psm_rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select(!`Spectrum number`)

  out$psm_sequence <- out$psm_sequence |> lapply(str_clean) |> unlist()
  out <- out |>
    dplyr::mutate(psm_peptide = purrr::map2(psm_sequence, Mods, omssa_peptide) |> unlist())

  return(out)
}

#' helper function to read in platform specific results
#'
#' @param sequence string
#' @param modifications string
#'
#' @return string
#'
omssa_peptide <- function(
    sequence = NULL,
    modifications = NULL
){

  sequence <- sequence |> stringr::str_extract_all('[A-Z]') |> unlist()
  modifications <- modifications |> stringr::str_split(",") |> unlist()

  masses <- modifications |> stringr::str_remove("\\:.+") |> lapply(omssa_mod) |> unlist()
  local <- modifications |> stringr::str_remove(".+\\:")

  out <- ''
  for(i in 1:length(sequence)){
    if(i %in% local){
      sequence[i] <- paste0("[", sequence[i], num_trunc(masses[which(i == local)],2), "]")
    }
    out <- paste0(out, sequence[i])
  }
  return(out)
}

omssa_mod <- function(x){
  switch (x,
          'Methylation of K' = {14.01565},
          'Oxidation of M' = {15.99491},
          'Carboxymethyl C' = {58.00547},
          'Carbamidomethyl C' = {57.02146},
          'Deamidation of N and Q' = {0.98401},
          'Propionamide C' = {71.03711},
          'Phosphorylation of S' = {79.96633},
          'Phosphorylation of T' = {79.96633},
          'Phosphorylation of Y' = {79.96633},
          'M cleavage from protein n-term' = {-131.04048},
          'Acetylation of protein n-term' = {42.01056},
          'Methylation of protein n-term' = {14.01565},
          'Tri-methylation of protein n-term' = {42.04695},
          'Beta methythiolation of D' = {45.98772},
          'Methylation of Q' = {14.01565},
          'Tri-methylation of K' = {42.04695},
          'Methylation of D' = {14.01565},
          'Methylation of E' = {14.01565},
          'Methylation of peptide c-term' = {14.01565},
          'Tri-deuteromethylation of D' = {17.03448},
          'Tri-deuteromethylation of E' = {17.03448},
          'Tri-deuteromethylation of peptide c-term' = {17.03448},
          'N-formyl met addition' = {159.03539},
          '2-amino-3-oxo-butanoic acid T' = {-2.01565},
          'Acetylation of K' = {42.01056},
          'Amidation of peptide c-term' = {-0.98401},
          'Beta-methylthiolation of D (duplicate of 13)' = {45.98772},
          'Carboxyamidomethylation of K' = {57.02146},
          'Carboxyamidomethylation of H' = {57.02146},
          'Carboxyamidomethylation of D' = {57.02146},
          'Carboxyamidomethylation of E' = {57.02146},
          'Carbamylation of K' = {43.00581},
          'Carbamylation of n-term peptide' = {43.00581},
          'Citrullination of R' = {0.98401},
          'Oxidation of C to cysteic acid' = {47.98474},
          'Di-iodination of Y' = {251.79329},
          'Di-methylation of K' = {28.03130},
          'Di-methylation of R' = {28.03130},
          'Di-methylation of peptide n-term' = {28.03130},
          'Oxidation of F to dihydroxyphenylalanine' = {31.98982},
          'Gammathiopropionylation of K' = {87.99828},
          'Gammathiopropionylation of peptide n-term' = {87.99828},
          'Farnesylation of C' = {204.18780},
          'Formylation of K' = {27.99491},
          'Formylation of peptide n-term' = {27.99491},
          'Oxidation of W to formylkynurenin' = {31.98982},
          'Fluorophenylalanine' = {17.99057},
          'Beta-carboxylation of D' = {43.98982},
          'Gamma-carboxylation of E' = {43.98982},
          'Geranyl-geranyl' = {272.25040},
          'Glucuronylation of protein n-term' = {176.03208},
          'Glutathione disulfide' = {305.06815},
          'ubiquitinylation residue' = {114.04292},
          'Guanidination of K' = {42.02179},
          'Oxidation of H to N' = {-23.01598},
          'Oxidation of H to D' = {-22.03196},
          'Homoserine' = {-29.99280},
          'Homoserine lactone' = {-48.00337},
          'Oxidation of W to hydroxykynurenin' = {19.98982},
          'Hydroxylation of D' = {15.99491},
          'Hydroxylation of K' = {15.99491},
          'Hydroxylation of N' = {15.99491},
          'Hydroxylation of P' = {15.99491},
          'Hydroxylation of F' = {15.99491},
          'Hydroxylation of Y' = {15.99491},
          'Iodination of Y' = {125.89664},
          'Oxidation of W to kynurenin' = {3.99491},
          'Lipoyl K' = {188.03295},
          'Methyl ester of peptide c-term (duplicate of 18)' = {14.01565},
          'Methyl ester of D' = {14.01565},
          'Methyl ester of E (duplicate of 17)' = {14.01565},
          'Methyl ester of S' = {14.01565},
          'Methyl ester of Y' = {14.01565},
          'Methyl C' = {14.01565},
          'Methyl H' = {14.01565},
          'Methyl N' = {14.01565},
          'Methylation of peptide n-term' = {14.01565},
          'Methyl R' = {14.01565},
          'Myristoleylation of G' = {208.18271},
          'Myristoyl-4H of G' = {206.16706},
          'Myristoylation of peptide n-term G' = {210.19836},
          'Myristoylation of K' = {210.19836},
          'Formylation of protein n-term' = {27.99491},
          'NEM C' = {125.04767},
          'NIPCAM' = {99.06841},
          'Oxidation of W to nitro' = {44.98507},
          'Oxidation of Y to nitro' = {44.98507},
          'O18 on peptide n-term' = {2.00424},
          'Di-O18 on peptide n-term' = {4.0084},
          'Oxidation of H' = {15.99491},
          'Oxidation of W' = {15.99491},
          'Phosphopantetheine S' = {340.08579},
          'Palmitoylation of C' = {238.22966},
          'Palmitoylation of K' = {238.22966},
          'Palmitoylation of S' = {238.22966},
          'Palmitoylation of T' = {238.22966},
          'Phosphorylation of S with prompt loss' = {-18.01056},
          'Phosphorylation of T with prompt loss' = {-18.01056},
          'Phosphorylation with prompt loss on Y' = {-18.01056},
          'Phosphorylation with neutral loss on C' = {79.96633},
          'Phosphorylation with neutral loss on D' = {79.96633},
          'Phosphorylation with neutral loss on H' = {79.96633},
          'Propionyl light K' = {56.02621},
          'Propionyl light on peptide n-term' = {56.02621},
          'Propionyl heavy K' = {59.03627},
          'Propionyl heavy peptide n-term' = {59.03627},
          'Pyridyl K' = {119.03711},
          'Pyridyl peptide n-term' = {119.03711},
          'Pyro-cmC' = {-17.02654},
          'Pyro-glu from n-term E' = {-18.01056},
          'Pyro-glu from n-term Q' = {-17.02654},
          'Oxidation of P to pyroglutamic acid' = {13.97926},
          'S-pyridylethylation of C' = {105.05784},
          'SeMet' = {47.94444},
          'Sulfation of Y' = {79.95681},
          'Sulphone of M' = {31.98982},
          'Tri-iodination of Y' = {377.68994},
          'Tri-methylation of R' = {42.04695},
          'N-acyl diglyceride cysteine' = {788.72577},
          'ICAT light' = {227.12699},
          'ICAT heavy' = {236.15718},
          'CAMthiopropanoyl K' = {145.01974},
          'Phosphorylation with neutral loss on S' = {79.96633},
          'Phosphorylation with neutral loss on T' = {79.96633},
          'Phosphorylation of S with ETD loss' = {79.96633},
          'Phosphorylation of T with ETD loss' = {79.96633},
          'Heavy arginine-13C6' = {6.02012},
          'Heavy arginine-13C6-15N4' = {10.00826},
          'Heavy lysine-13C6' = {6.02012},
          'PNGasF in O18 water' = {2.98826},
          'Beta elimination of S' = {-18.01056},
          'Beta elimination of T' = {-18.01056},
          'Oxidation of C to sulfinic acid' = {31.98982},
          'Arginine to ornithine' = {-42.02179},
          'Dehydro of S and T' = {-18.01056},
          'Carboxykynurenin of W' = {47.9847438},
          'Sumoylation of K' = {484.228},
          'ITRAQ114 on nterm' = {144.10591},
          'ITRAQ114 on K' = {144.10591},
          'ITRAQ114 on Y' = {144.10591},
          'ITRAQ115 on nterm' = {144.09959},
          'ITRAQ115 on K' = {144.09959},
          'ITRAQ115 on Y' = {144.09959},
          'ITRAQ116 on nterm' = {144.10206},
          'ITRAQ116 on K' = {144.10206},
          'ITRAQ116 on Y' = {144.10206},
          'ITRAQ117 on nterm' = {144.10206},
          'ITRAQ117 on K' = {144.10206},
          'ITRAQ117 on Y' = {144.10206},
          'MMTS on C' = {45.98772},
          'Heavy lysine - 2H4' = {4.02510},
          'Heavy lysine - 13C6 15N2' = {8.01419},
          'Asparagine HexNAc' = {203.07937},
          'Asparagine dHexHexNAc' = {349.13728},
          'Serine HexNAc' = {203.07937},
          'Threonine HexNAc' = {203.07937},
          'Palmitoleyl of S' = {236.21401},
          'Palmitoleyl of C' = {236.21401},
          'Palmitoleyl of T' = {236.21401},
          'CHD2-di-methylation of K' = {32.05640},
          'CHD2-di-methylation of peptide n-term' = {32.05640},
          'Maleimide-PEO2-Biotin of C' = {525.22571},
          'Phosphorylation of H' = {79.96633},
          'Oxidation of C' = {15.99491},
          'Oxidation of Y (duplicate of 64)' = {15.99491},
          'Uniblue A on K' = {484.03989},
          'Deamidation of N' = {0.98401},
          'Trideuteration of L (SILAC)' = {3.01883},
          'TMT duplex on K (old)' = {225.15583},
          'TMT duplex on n-term peptide (old)' = {225.15583},
          'TMT 6-plex on K (old)' = {229.16293},
          'TMT 6-plex on n-term peptide (old)' = {229.16293},
          'ITRAQ8plex:13C(7)15N(1) on nterm' = {304.20536},
          'ITRAQ8plex:13C(7)15N(1) on K' = {304.20536},
          'ITRAQ8plex:13C(7)15N(1) on Y' = {304.20536},
          'ITRAQ8plex:13C(6)15N(2) on nterm' = {304.19904},
          'ITRAQ8plex:13C(6)15N(2) on K' = {304.19904},
          'ITRAQ8plex:13C(6)15N(2) on Y' = {304.19904},
          'Selenocysteine' = {47.94444},
          'Carboxymethylated selenocysteine' = {105.94992},
          {
            0.0
          })

}
