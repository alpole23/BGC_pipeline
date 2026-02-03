process BIGSLICE {
    tag "$taxon"
    label 'process_high'
    publishDir "${params.outdir}/bigslice_results", mode: 'copy'
    cache 'lenient'  // Use lenient caching for directory inputs

    input:
    val taxon
    path antismash_results
    path bigslice_models
    path taxonomy_map
    path name_map

    output:
    path "${Utils.sanitizeTaxon(taxon)}/", emit: bigslice_dir

    script:
    def threshold = params.bigslice_threshold ?: "900"
    def safe_taxon = Utils.sanitizeTaxon(taxon)
    def query_mode = params.bigslice_query_mode ?: false
    def query_db = params.bigslice_query_db ?: ""
    def query_flag = query_mode && query_db ? "--query ${query_db}" : ""
    """
    echo "Running BiG-SLiCE on ${taxon}"
    echo "Using ${task.cpus} cores"
    echo "Query mode: ${query_mode}"

    # Test BiG-SLiCE installation
    bigslice --version

    # Set BiG-SLiCE models path
    export MODELS_PATH=\$(readlink -f ${bigslice_models})
    echo "Using BiG-SLiCE models from: \$MODELS_PATH"

    # Create output directory
    mkdir -p ${safe_taxon}

    # Create input directory with proper nested structure for BiG-SLiCE
    # Structure: antismash_input/dataset_name/genome_folder/gbk_files
    DATASET_NAME="${safe_taxon}"
    mkdir -p antismash_input/\$DATASET_NAME
    mkdir -p antismash_input/taxonomy

    # Process each antiSMASH result directory and organize by genome
    GENOME_COUNT=0
    for dir in ${antismash_results}; do
        # Resolve symlinks to actual directory path
        REAL_DIR=\$(readlink -f "\$dir")
        if [ -d "\$REAL_DIR" ]; then
            # Extract genome name from the directory path (last component)
            GENOME_NAME=\$(basename "\$REAL_DIR")
            echo "Processing genome: \$GENOME_NAME"

            # Create genome-specific folder within dataset
            GENOME_DIR="antismash_input/\$DATASET_NAME/\$GENOME_NAME"
            mkdir -p "\$GENOME_DIR"

            # First try to find region-specific GenBank files
            REGION_FILES=\$(find "\$REAL_DIR" -name "*.region*.gbk" 2>/dev/null | wc -l)

            if [ "\$REGION_FILES" -gt 0 ]; then
                echo "  Found \$REGION_FILES region GenBank files"
                find "\$REAL_DIR" -name "*.region*.gbk" -exec ln -s {} "\$GENOME_DIR/" \\;
                GENOME_COUNT=\$((GENOME_COUNT + 1))
            else
                # If no region files, use the complete GenBank file
                MAIN_GBK=\$(find "\$REAL_DIR" -maxdepth 1 -name "*.gbk" -type f | head -1)
                if [ -n "\$MAIN_GBK" ]; then
                    echo "  Using complete GenBank file: \$(basename \$MAIN_GBK)"
                    ln -s "\$MAIN_GBK" "\$GENOME_DIR/"
                    GENOME_COUNT=\$((GENOME_COUNT + 1))
                fi
            fi
        fi
    done

    echo "Total genomes organized: \$GENOME_COUNT"

    if [ "\$GENOME_COUNT" -eq 0 ]; then
        echo "ERROR: No GenBank files found"
        echo "Listing what we received:"
        ls -la
        exit 1
    fi

    # Generate taxonomy file from taxonomy_map.json
    echo "Generating taxonomy file..."
    python3 << PYEOF
import json
import os

# Read taxonomy_map.json
with open('${taxonomy_map}', 'r') as f:
    taxonomy_data = json.load(f)

# Read name_map.json to get renamed genome names
with open('${name_map}', 'r') as f:
    name_map_data = json.load(f)

# BiG-SLiCE taxonomy format (tab-separated):
# Genome folder	Kingdom	Phylum	Class	Order	Family	Genus	Species	Organism
taxonomy_file = 'antismash_input/taxonomy/\${DATASET_NAME}_taxonomy.tsv'

with open(taxonomy_file, 'w') as out:
    # Write header (required by BiG-SLiCE)
    out.write('# Genome folder\\tKingdom\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\tOrganism\\n')

    # Write taxonomy for each genome
    for assembly_id, data in taxonomy_data.items():
        lineage = data.get('lineage', {})
        organism = data.get('organism', '')

        # Get renamed genome name from name_map.json
        # This matches the actual folder name in antismash_input/dataset_name/
        if assembly_id in name_map_data:
            genome_folder = name_map_data[assembly_id]
        else:
            # Fallback to assembly_id if not found in name_map
            genome_folder = assembly_id.replace('.', '_')
            print(f"Warning: Assembly {assembly_id} not found in name_map, using {genome_folder}")

        # Get taxonomic ranks (BiG-SLiCE expects these specific ranks)
        kingdom = lineage.get('kingdom', {}).get('name', '')
        phylum = lineage.get('phylum', {}).get('name', '')
        tax_class = lineage.get('class', {}).get('name', '')
        order = lineage.get('order', {}).get('name', '')
        family = lineage.get('family', {}).get('name', '')
        genus = lineage.get('genus', {}).get('name', '')
        species = lineage.get('species', {}).get('name', '')

        # Write taxonomy line
        # Genome folder path should end with '/' as per BiG-SLiCE requirement
        out.write(f'{genome_folder}/\\t{kingdom}\\t{phylum}\\t{tax_class}\\t{order}\\t{family}\\t{genus}\\t{species}\\t{organism}\\n')

print(f"Generated taxonomy file: {taxonomy_file}")
print(f"Total genomes in taxonomy: {len(taxonomy_data)}")
PYEOF

    echo "Taxonomy file created. Verifying..."
    ls -la antismash_input/taxonomy/
    cat antismash_input/taxonomy/\${DATASET_NAME}_taxonomy.tsv | head -5

    # Create datasets.tsv metadata file required by BiG-SLiCE
    # Format: dataset_name<tab>path_to_folder<tab>taxonomy_path<tab>description
    # IMPORTANT: Must include header line starting with '#'
    echo "Creating datasets.tsv..."
    echo -e "# Dataset name\tPath to folder\tPath to taxonomy\tDescription" > antismash_input/datasets.tsv
    echo -e "\$DATASET_NAME\t\$DATASET_NAME/\ttaxonomy/\${DATASET_NAME}_taxonomy.tsv\t${taxon}" >> antismash_input/datasets.tsv

    echo "datasets.tsv contents:"
    cat antismash_input/datasets.tsv

    echo "Directory structure:"
    ls -lR antismash_input/ | head -50

    # Run BiG-SLiCE
    # BiG-SLiCE creates a database and performs clustering in a single run
    # If query_mode is enabled, compares against existing database instead of clustering
    bigslice \\
        -i antismash_input \\
        ${safe_taxon}/bigslice_output \\
        --threshold ${threshold} \\
        --program_db_folder \$MODELS_PATH \\
        -t ${task.cpus} \\
        ${query_flag}

    # Fix BiG-SLiCE database to support antiSMASH 8.x
    # BiG-SLiCE only includes type definitions for antiSMASH 4, 5, and MIBiG
    # We need to add 'as8' type to support antiSMASH 8.x GenBank files
    echo "Adding antiSMASH 8 support to BiG-SLiCE database..."
    DB_PATH="${safe_taxon}/bigslice_output/result/data.db"

    # Check if as8 type already exists (to make this idempotent)
    AS8_EXISTS=\$(sqlite3 "\$DB_PATH" "SELECT COUNT(*) FROM enum_bgc_type WHERE code='as8';")

    if [ "\$AS8_EXISTS" -eq 0 ]; then
        sqlite3 "\$DB_PATH" "INSERT INTO enum_bgc_type (code, description) VALUES ('as8', 'antiSMASH8 regionXXX.gbk');"
        echo "Added 'as8' type to enum_bgc_type table"
    else
        echo "'as8' type already exists in enum_bgc_type table"
    fi

    # Create empty reports database to prevent web interface errors
    # The reports feature is for query mode, which we don't use in standard clustering
    # Creating an empty database prevents the Reports page from crashing
    echo "Creating empty reports database structure..."
    REPORTS_DIR="${safe_taxon}/bigslice_output/reports"
    mkdir -p "\$REPORTS_DIR"
    REPORTS_DB="\$REPORTS_DIR/reports.db"

    if [ ! -f "\$REPORTS_DB" ]; then
        sqlite3 "\$REPORTS_DB" "CREATE TABLE IF NOT EXISTS reports (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            type TEXT,
            creation_date TEXT,
            name TEXT
        );"
        echo "Created empty reports database at \$REPORTS_DB"
    else
        echo "Reports database already exists"
    fi

    echo "BiG-SLiCE analysis complete"
    echo "Output directory: ${safe_taxon}/bigslice_output"
    """
}
