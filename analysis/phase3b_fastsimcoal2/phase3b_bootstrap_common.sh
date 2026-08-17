#!/bin/bash
################################################################################
# Shared helpers for fastsimcoal2 parametric bootstrap (Step 4)
################################################################################

extract_max_est_lhood() {
    local bestlhoods_file=$1
    awk '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                key = tolower($i)
                if (key == "maxestlhood" || key ~ /^maxestlhood/) {
                    idx = i
                }
            }
            next
        }
        NF > 0 { last = $0 }
        END {
            if (last == "") {
                exit 1
            }
            n = split(last, fields, /[[:space:]]+/)
            if (idx > 0 && idx <= n) {
                print fields[idx]
            } else if (n >= 2) {
                print fields[n - 1]
            } else {
                exit 1
            }
        }
    ' "${bestlhoods_file}"
}

# Convert a *_maxL.par (FREQ/OUTEXP estimation file) into a DNA simulation file.
# fastsimcoal2 2.8 cannot simulate SFS from FREQ 1 0 .par files with -d -s 0;
# use many independent DNA loci as in the official parametric bootstrap workflow.
prepare_bootstrap_sim_par() {
    local src_par=$1
    local dst_par=$2
    local n_loci=${3:-200000}
    local dna_length=${4:-100}

    awk -v n_loci="${n_loci}" -v dna_length="${dna_length}" '
        /^\/\/Number of independent loci/ {
            print
            getline
            print n_loci " 0"
            next
        }
        /^FREQ / {
            mut_rate = $4
            for (i = 5; i <= NF; i++) {
                extra = extra " " $i
            }
            print "DNA " dna_length " 0 " mut_rate extra
            next
        }
        { print }
    ' "${src_par}" > "${dst_par}"
}

find_simulated_sfs() {
    local rep_dir=$1
    local boot_prefix=$2

    local candidate="${rep_dir}/${boot_prefix}/${boot_prefix}_1/${boot_prefix}_DAFpop0.obs"
    if [ -f "${candidate}" ]; then
        echo "${candidate}"
        return 0
    fi

    candidate=$(find "${rep_dir}/${boot_prefix}" -name '*_DAFpop0.obs' -type f 2>/dev/null | head -1)
    if [ -n "${candidate}" ] && [ -f "${candidate}" ]; then
        echo "${candidate}"
        return 0
    fi

    return 1
}
