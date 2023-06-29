// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_rmstandem_extendr(void *dll);

void R_init_rmstandem(void *dll) {
    R_init_rmstandem_extendr(dll);
}
