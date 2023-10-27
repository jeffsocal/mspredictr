# parse test for SearchGUI xtandem results

load_all()
path <- "~/Local/data/project_jester/results/simulation_ecoli_10kPeptides/parse_test"
tbl <- glue::glue("{path}/peptides10000_simPTM_ecoli-BL21DE3-4156_1.t.xml") |>
  read_psms(platform = 'xtandem', cpus = 1)
