import os
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run_single_sample_mapping(sample, cpus, memory, time, project_input, project_output):
    """
    Maps reads to the single-sample assembly using bwa mem.
    """
    print(f"Running Single-Sample Mapping for sample: {sample}")

    # Define output directories (Separate for each sample)
    output_dir = f"{project_output}/Mapping/Single_Sample/{sample}"
    ensure_directory_exists(output_dir)

    # Reference and Input Files
    assembly_dir = f"{project_output}/Assembly/IDBA"
    contigs = f"{assembly_dir}/{sample}_idba/scaffold.fa"
    reads_dir = f"{project_output}/Host_Removal"
    r1 = f"{reads_dir}/{sample}_unmapped_R1.fastq.gz"
    r2 = f"{reads_dir}/{sample}_unmapped_R2.fastq.gz"
    mapped_bam = f"{output_dir}/{sample}_mapped.bam"
    sorted_bam = f"{output_dir}/{sample}_sorted.bam"
    indexed_bam = f"{sorted_bam}.bai"
    depth_file = f"{output_dir}/depth_{sample}.txt"

    jgi_summarize = get_tool_path("jgi_summarize")
    if not jgi_summarize:
        print(f"[ERROR] Jgi_Summarize tool not found in configuration")
        return
    bwa = get_tool_path("bwa")
    if not bwa:
        print(f"[ERROR] Bwa tool not found in configuration")
        return
    samtools = get_tool_path("samtools")
    if not samtools:
        print(f"[ERROR] Samtools tool not found in configuration")
        return

    # **Step 1: Index the reference scaffold if missing**
    index_files = [f"{contigs}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    if not all(os.path.exists(f) for f in index_files):
        print(f"[INFO] Indexing reference scaffold for sample {sample}")
        run_command(f"{bwa} index {contigs}")
    else:
        print(f"[INFO] Reference index files already exist for {sample}, skipping indexing.")

    # **Step 2: Run mapping only if BAM does not exist**
    if not os.path.exists(mapped_bam):
        print(f"[INFO] Mapping reads for sample {sample}")
        map_command = f"{bwa} mem -t {cpus} {contigs} {r1} {r2} | {samtools} fixmate -m -O bam - {mapped_bam}"
        run_command(map_command)
    else:
        print(f"[INFO] Mapped BAM file already exists for {sample}, skipping mapping.")

    # **Step 3: Sort BAM if not sorted**
    if not os.path.exists(sorted_bam):
        print(f"[INFO] Sorting BAM file for sample {sample}")
        sort_command = f"{samtools} sort {mapped_bam} -o {sorted_bam}"
        run_command(sort_command)
    else:
        print(f"[INFO] Sorted BAM file already exists for {sample}, skipping sorting.")

    # **Step 4: Index BAM if not already indexed**
    if not os.path.exists(indexed_bam):
        print(f"[INFO] Indexing BAM file for sample {sample}")
        index_command = f"{samtools} index {sorted_bam}"
        run_command(index_command)
    else:
        print(f"[INFO] BAM index already exists for {sample}, skipping indexing.")

    # **Step 5: Run depth summarization only if output file is missing**
    if not os.path.exists(depth_file):
        print(f"[INFO] Running depth summarization for sample {sample}")
        depth_command = f"bash -c 'source /usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate metabat_env && {jgi_summarize} --outputDepth {depth_file} {sorted_bam}'"
        run_command(depth_command)
    else:
        print(f"[INFO] Depth summarization already completed for {sample}, skipping.")


def run_coassembly_mapping(samples, cpus, memory, time, project_input, project_output ):
    """
    Maps reads to the co-assembly using bwa mem.
    """
    print(f"Running Co-Assembly Mapping for samples: {', '.join(samples)}")

    # Define output directories (Separate for each sample)
    for sample in samples:
        sample_output_dir = f"{project_output}/Mapping/Coassembly/{sample}"
        ensure_directory_exists(sample_output_dir)

    coassembly_dir = f"{project_output}/Assembly/MEGAHIT"
    contigs = f"{coassembly_dir}/coassembly_contigs.fa"
    reads_dir = f"{project_output}/Host_Removal"
    depth_file = f"{project_output}/Mapping/Coassembly/depth_all.txt"

    jgi_summarize = get_tool_path("jgi_summarize")
    bwa = get_tool_path("bwa")
    samtools = get_tool_path("samtools")

    mapped_bams = []

    # **Step 1: Index the co-assembly reference if missing**
    index_files = [f"{contigs}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    if not all(os.path.exists(f) for f in index_files):
        print(f"[INFO] Indexing co-assembly reference")
        run_command(f"{bwa} index {contigs}")
    else:
        print(f"[INFO] Co-assembly reference index already exists, skipping indexing.")

    for sample in samples:
        r1 = f"{reads_dir}/{sample}_unmapped_R1.fastq.gz"
        r2 = f"{reads_dir}/{sample}_unmapped_R2.fastq.gz"
        sample_output_dir = f"{project_output}/Mapping/Coassembly/{sample}"
        mapped_bam = f"{sample_output_dir}/{sample}_coassembly_mapped.bam"
        sorted_bam = f"{sample_output_dir}/{sample}_coassembly_sorted.bam"
        mapped_bams.append(sorted_bam)

        # **Step 2: Run mapping only if BAM does not exist**
        if not os.path.exists(mapped_bam):
            print(f"[INFO] Mapping reads for sample {sample}")
            map_command = f"{bwa} mem -t {cpus} {contigs} {r1} {r2} | {samtools} fixmate -m -O bam - {mapped_bam}"
            run_command(map_command)
        else:
            print(f"[INFO] Mapped BAM already exists for {sample}, skipping mapping.")

        # **Step 3: Sort BAM if not sorted**
        if not os.path.exists(sorted_bam):
            print(f"[INFO] Sorting BAM file for sample {sample}")
            sort_command = f"{samtools} sort {mapped_bam} -o {sorted_bam}"
            run_command(sort_command)
        else:
            print(f"[INFO] Sorted BAM file already exists for {sample}, skipping sorting.")

        # **Step 4: Index BAM if not already indexed**
        if not os.path.exists(f"{sorted_bam}.bai"):
            print(f"[INFO] Indexing BAM file for sample {sample}")
            index_command = f"{samtools} index {sorted_bam}"
            run_command(index_command)
        else:
            print(f"[INFO] BAM index already exists for {sample}, skipping indexing.")

    # **Step 5: Run depth summarization only if missing**
    if not os.path.exists(depth_file):
        print(f"[INFO] Running depth summarization for co-assembly")
        depth_command = f"bash -c 'source /usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate metabat_env && {jgi_summarize} --outputDepth {depth_file} {' '.join(mapped_bams)}'"
        run_command(depth_command)
    else:
        print(f"[INFO] Depth summarization for co-assembly already exists, skipping.")
