# parse test for SearchGUI xtandem results

load_all()
path <- "~/Local/data/project_jester/results/simulation_ecoli_10kPeptides/parse_test"

tbl <- list(
  glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.sage.tsv") |>
    read_psms(platform = 'sage', cpus = 1),

  glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.ms-amanda.csv") |>
    read_psms(platform = 'ms_amanda', cpus = 1),

  glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.comet.txt") |>
    read_psms(platform = 'comet', cpus = 1),

  glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.csv") |>
    read_psms(platform = 'omssa', cpus = 1),

  glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.t.xml") |>
    read_psms(platform = 'xtandem', cpus = 1)
)
