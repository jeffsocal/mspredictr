# parse test for SearchGUI comet results

load_all()
path <- "~/Local/data/project_jester/results/simulation_ecoli_10kPeptides/parse_test"
tbl <- glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.csv") |>
  read_psms(platform = 'omssa', cpus = 1)
