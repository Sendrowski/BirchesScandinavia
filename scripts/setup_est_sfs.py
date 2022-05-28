from snakemake.shell import shell

try:
    out = snakemake.output[0]
    tmp_dir = snakemake.resources.tmpdir
    max_sites = snakemake.config['est_sfs_max_sites']
except NameError:
    # testing
    out = "scratch/est-sfs"
    tmp_dir = "/tmp"
    max_sites = 100000

# set up est-sfs
# Note: the compilation only works on Linux
shell(f"""
    set -x
    work_dir=$(realpath .)
    cd {tmp_dir}
    
    wget https://sourceforge.net/projects/est-usfs/files/est-sfs-release-2.03.tar.gz/download -O est-sfs.tar.gz
    
    rm -rf est-sfs
    mkdir est-sfs
    tar -xvf est-sfs.tar.gz -C est-sfs --strip-components 1
    
    cd est-sfs
    sed -i 's/#define max_config 100000/#define max_config {max_sites}/g' est-sfs.c
    export C_INCLUDE_PATH=$CONDA_PREFIX/include
    make
    
    mv est-sfs "$work_dir/{out}"
""")
