process COLLECT_VERSIONS {
    tag "versions"
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path antismash_done  // Dependency to ensure conda envs exist
    path bigscape_done   // Optional dependency for BiG-SCAPE
    path gtdbtk_done     // Optional dependency for GTDB-Tk

    output:
    path "software_versions.json", emit: versions

    script:
    def nf_version = workflow.nextflow.version
    def work_conda = "${workDir}/conda"
    """
    #!/usr/bin/env python3
    import json
    import subprocess
    import os
    import re
    from pathlib import Path

    def find_executable(name, search_dir):
        '''Find executable in conda environments.'''
        for env_dir in Path(search_dir).glob("env-*"):
            bin_path = env_dir / "bin" / name
            if bin_path.exists():
                return str(bin_path)
        return None

    def get_version(executable, pattern=None, version_arg="--version"):
        '''Run version command and extract version.'''
        try:
            result = subprocess.run(
                [executable, version_arg],
                capture_output=True, text=True, timeout=30
            )
            output = result.stdout + result.stderr
            if pattern:
                match = re.search(pattern, output)
                if match:
                    return match.group(1)
            # Return first non-empty line
            for line in output.strip().split("\\n"):
                if line.strip():
                    return line.strip()
            return "unknown"
        except Exception as e:
            return f"error: {str(e)}"

    work_conda = "${work_conda}"
    versions = {}

    # antiSMASH
    antismash_bin = find_executable("antismash", work_conda)
    if antismash_bin:
        versions["antismash"] = get_version(antismash_bin, r"antiSMASH\\s+([\\d.]+)")
    else:
        versions["antismash"] = "not installed"

    # BiG-SCAPE - check if it was run (bigscape_done is not a placeholder)
    bigscape_ran = os.path.isdir("${bigscape_done}") or (os.path.isfile("${bigscape_done}") and "NO_BIGSCAPE" not in "${bigscape_done}")
    bigscape_bin = find_executable("bigscape", work_conda)
    if bigscape_bin:
        versions["bigscape"] = get_version(bigscape_bin, r"BiG-SCAPE\\s+([\\d.]+)")
    elif bigscape_ran:
        versions["bigscape"] = "used (version unknown)"
    else:
        versions["bigscape"] = "not used"

    # BiG-SLiCE - check named env first, then work conda
    bigslice_version = "not used"
    try:
        result = subprocess.run(
            ["conda", "run", "-n", "bigslice", "bigslice", "--version"],
            capture_output=True, text=True, timeout=30
        )
        match = re.search(r"version\\s+([\\d.]+)", result.stdout + result.stderr)
        if match:
            bigslice_version = match.group(1)
    except:
        bigslice_bin = find_executable("bigslice", work_conda)
        if bigslice_bin:
            bigslice_version = get_version(bigslice_bin, r"version\\s+([\\d.]+)")
    versions["bigslice"] = bigslice_version

    # GTDB-Tk - check if it was run (gtdbtk_done is not a placeholder)
    gtdbtk_ran = os.path.isdir("${gtdbtk_done}") or (os.path.isfile("${gtdbtk_done}") and "NO_GTDBTK" not in "${gtdbtk_done}")
    gtdbtk_bin = find_executable("gtdbtk", work_conda)
    if gtdbtk_bin:
        versions["gtdbtk"] = get_version(gtdbtk_bin, r"gtdbtk:\\s*v?([\\d.]+)")
    elif gtdbtk_ran:
        versions["gtdbtk"] = "used (version unknown)"
    else:
        versions["gtdbtk"] = "not used"

    # TaxonKit (uses 'version' subcommand, not '--version' flag)
    taxonkit_bin = find_executable("taxonkit", work_conda)
    if taxonkit_bin:
        versions["taxonkit"] = get_version(taxonkit_bin, r"taxonkit\\s+v?([\\d.]+)", "version")
    else:
        versions["taxonkit"] = "not installed"

    # Nextflow and pipeline
    versions["nextflow"] = "${nf_version}"
    versions["pipeline_version"] = "1.0.0"

    with open("software_versions.json", "w") as f:
        json.dump(versions, f, indent=2)

    print("Software versions collected:")
    for tool, version in versions.items():
        print(f"  {tool}: {version}")
    """
}
