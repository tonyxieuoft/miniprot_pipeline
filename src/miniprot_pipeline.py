import os
import Bio.Align

N_to_AA = {
"TTT": "F",
"TTC": "F",
"TTA": "L",
"TTG": "L",
"TCT": "S",
"TCC": "S",
"TCA": "S",
"TCG": "S",
"TAT": "Y",
"TAC": "Y",
"TAA": "*",
"TAG": "*",
"TGT": "C",
"TGC": "C",
"TGA": "*",
"TGG": "W",
"CTT": "L",
"CTC": "L",
"CTA": "L",
"CTG": "L",
"CCT": "P",
"CCC": "P",
"CCA": "P",
"CCG": "P",
"CAT": "H",
"CAC": "H",
"CAA": "Q",
"CAG": "Q",
"CGT": "R",
"CGC": "R",
"CGA": "R",
"CGG": "R",
"ATT": "I",
"ATC": "I",
"ATA": "I",
"ATG": "M",
"ACT": "T",
"ACC": "T",
"ACA": "T",
"ACG": "T",
"AAT": "N",
"AAC": "N",
"AAA": "K",
"AAG": "K",
"AGT": "S",
"AGC": "S",
"AGA": "R",
"AGG": "R",
"GTT": "V",
"GTC": "V",
"GTA": "V",
"GTG": "V",
"GCT": "A",
"GCC": "A",
"GCA": "A",
"GCG": "A",
"GAT": "D",
"GAC": "D",
"GAA": "E",
"GAG": "E",
"GGT": "G",
"GGC": "G",
"GGA": "G",
"GGG": "G"}

def listdir(dir_path):
    """
    Lists the names of files contained within the directory specified in dir_path. Distinct from "os.listdir" in that
    it skips over files beginning with ".".
    :param dir_path: Path to a directory.
    :return: Names of files in the directory.
    """

    paths_to_return = []

    for f in os.listdir(dir_path):
        if not f.startswith("."):
            paths_to_return.append(f)

    return paths_to_return


def make_unique_directory(path):
    """
    Creates a new directory specified by path. If such a directory already exists, appends a version number to the end.

    :param path: Path to create the new directory
    :return: Path of the created directory. If a directory already exists at the original path, returns the new path
    with version number appended.
    """

    if os.path.isfile(path):

        counter = 1
        while os.path.isfile(path + str(counter)):
            counter += 1

        path = path + str(counter)

    os.mkdir(path)
    return path


def reference_sequences_nuc_to_protein(reference_sequences_nuc, download_dir):

    reference_sequences_prot = make_unique_directory(os.path.join(download_dir, "reference_sequences_protein"))

    for taxa_dir in listdir(reference_sequences_nuc):

        taxa_nuc_path = os.path.join(reference_sequences_nuc, taxa_dir)
        taxa_prot_path = make_unique_directory(os.path.join(reference_sequences_prot, taxa_dir))

        for gene_nuc_file in listdir(taxa_nuc_path):

            gene_nuc_path = os.path.join(taxa_nuc_path, gene_nuc_file)
            f = open(gene_nuc_path, "r")

            N_sequence = ""

            first_line = f.readline()

            line = first_line
            while line != "":
                if len(line.strip()) > 0:
                    if line.strip()[0] != ">":
                        N_sequence += line.strip()

                line = f.readline()

            AA_sequence = ""
            if len(N_sequence) % 3 != 0:
                print("there's a problem with " + gene_nuc_file)
            else:
                for i in range(len(N_sequence) // 3):
                    AA_sequence += N_to_AA[N_sequence[3 * i:3 * i + 3]]

            gene_prot_path = os.path.join(taxa_prot_path, gene_nuc_file)
            f2 = open(gene_prot_path, "w")
            f2.write(first_line + AA_sequence + "\n")

    return reference_sequences_prot


def miniprot_index(master_genome_fastas, download_path):

    """
    Create miniprot genome indexes for each genome fasta file in master_genome_fastas.

    :param master_genome_fastas: Directory containing genome fasta files. Must be formatted in a particular manner:
    the master directory must contain separate child directories for each taxon of interest, and the child directories
    must contain only fasta files named by species.
    :param download_path: Path to download the newly created miniprot genome indexes.
    :return: Path to the master directory containing the newly created miniprot genome indexes. Output folder
    structure mirrors the input folder structure.
    """

    master_genome_index_path = make_unique_directory(os.path.join(download_path, "master_genome_indexes"))

    for taxa_genome_fastas in listdir(master_genome_fastas):

        taxa_genome_index_dir_path = os.path.join(master_genome_index_path, taxa_genome_fastas)
        os.mkdir(taxa_genome_index_dir_path)

        taxa_genome_fastas_dir_path = os.path.join(master_genome_fastas, taxa_genome_fastas)
        for genome_fasta in listdir(taxa_genome_fastas_dir_path):

            species_name = os.path.splitext(genome_fasta)[0]

            command = ("nice -3 /crun2/storage5/AnthonyR/miniprot/miniprot" +
                       " -t 8" + " -d " +
                       os.path.join(taxa_genome_index_dir_path, species_name + ".mpi") + " " +  # not sure how to make sure it produces format "species_name.mpi" here
                       os.path.join(download_path, taxa_genome_fastas_dir_path, genome_fasta))

            print(command)
            os.system(command)

    return master_genome_index_path

            # miniprot aligning - aligns reference protein fastas to genome and produces gff file with cds info

def miniprot_align(reference_sequences, master_genome_indexes, download_path):  # again, can remove the out path if needed

    """
    Make sequence homology-based gene predictions through miniprot.

    :param reference_sequences: Directory containing (protein-level) reference sequences to base predictions off of.
    The master directory must contain child directories separating each taxon of interest. Each child directory must
    contain only .faa files named by gene.
    :param master_genome_indexes: Directory containing genome indexes created through the "miniprot_index" function
    :param download_path: Path to download the gene predictions.
    :return: Path to output directory containing predicted sequences in .gff format, which specifies ranges of coding
    sequences of the gene(s) within the genome assembly. Format of the output directory matches the format of the
    reference sequence input directory.
    """

    # make a master directory for miniprot alignments
    master_alignment_dir_path = make_unique_directory(os.path.join(download_path, "master_miniprot_alignment"))

    print(reference_sequences)
    print(listdir(reference_sequences))

    # taxa directory level of reference sequence directory
    for reference_taxa_dir in listdir(reference_sequences):
        reference_taxa_dir_path = os.path.join(reference_sequences, reference_taxa_dir)

        # make taxa directory level for miniprot alignments
        alignment_taxa_dir_path = os.path.join(master_alignment_dir_path, reference_taxa_dir)
        os.mkdir(alignment_taxa_dir_path)

        # go through the AAs contained in each file within the taxa directory level of ref. seq directory
        for protein_faa in listdir(reference_taxa_dir_path):

            # get the name of the AA --> needed for naming the alignment folder (level 3)
            protein_faa_path = os.path.join(reference_taxa_dir_path, protein_faa)
            protein_name = os.path.splitext(protein_faa)[0]

            # make gene directory level for miniprot alignments
            alignment_gene_dir_path = os.path.join(alignment_taxa_dir_path, protein_name)
            os.mkdir(alignment_gene_dir_path)

            # this can be done since the taxa dirs are named the name across the master directories
            genome_index_taxa_path = os.path.join(master_genome_indexes, reference_taxa_dir)
            for genome_index in listdir(genome_index_taxa_path):

                species_name = os.path.splitext(genome_index)[0]

                command = ("nice -3 /crun2/storage5/AnthonyR/miniprot/miniprot"
                          " -t 8 -I" +
                          " --gff " + os.path.join(genome_index_taxa_path, genome_index) + " " +
                          protein_faa_path +
                          " > " + os.path.join(alignment_gene_dir_path, species_name + ".gff"))

                print(command)
                os.system(command)

    return master_alignment_dir_path

                  # agat gff cds extraction - uses gff file from last step to extract cds from genome assembly

def agat_gff_cds_extraction(master_genome_fastas, master_miniprot_alignments, download_path):

    """
    Extracts nucleotide sequences from genome assemblies for the genes predicted by miniprot.

    :param master_genome_fastas: Directory containing genome assemblies in .fasta format. These genome assemblies must
    be the same assemblies that miniprot based its predictions off of.
    :param master_miniprot_alignments: Directory containing miniprot gene predictions in .gff format.
    :param download_path: Path to download the extracted nucleotide sequences.
    :return: Path to the output directory containing the extracted nucleotide sequences. The format of the output
    directory mirrors the structure of the input .gff file-containing directory.
    """


    master_agat_output_dir_path = make_unique_directory(os.path.join(download_path, "master_agat_output_dir_updated"))

    for miniprot_alignment_taxa_dir in listdir(master_miniprot_alignments):
        miniprot_alignment_taxa_dir_path = os.path.join(master_miniprot_alignments, miniprot_alignment_taxa_dir)

        agat_taxa_dir_path = os.path.join(master_agat_output_dir_path, miniprot_alignment_taxa_dir)
        os.mkdir(agat_taxa_dir_path)

        for miniprot_alignment_gene_dir in listdir(miniprot_alignment_taxa_dir_path):
            miniprot_alignment_gene_dir_path = os.path.join(miniprot_alignment_taxa_dir_path, miniprot_alignment_gene_dir)

            agat_gene_dir_path = os.path.join(agat_taxa_dir_path, miniprot_alignment_gene_dir)
            os.mkdir(agat_gene_dir_path)

            for gff_file in listdir(miniprot_alignment_gene_dir_path):

                gff_file_path = os.path.join(miniprot_alignment_gene_dir_path, gff_file)

                species_name = os.path.splitext(gff_file)[0]
                genome_fasta_path =  os.path.join(master_genome_fastas, miniprot_alignment_taxa_dir, species_name + ".fna")

                command = ("nice -3 agat_sp_extract_sequences.pl"
                          " -g " + gff_file_path +
                          " -f " + genome_fasta_path +
                          " -o " + os.path.join(agat_gene_dir_path, os.path.splitext(gff_file)[0] + ".fna"))

                print(command)
                os.system(command)

    return master_agat_output_dir_path

def combine_agat_files_by_gene(master_agat_output_dir, download_dir):

    """
    Combine all miniprot predictions for a single gene into one .fasta file.

    :param master_agat_output_dir: Directory containing gene predictions in nucleotide sequence format
    :param download_dir: Directory to download combined predictions.
    :return: Path to directory containing .fasta files with predictions combined by gene.
    """

    combined_agat_file_dir = make_unique_directory(os.path.join(download_dir, "combined_agat_files"))

    for taxa_dir in listdir(master_agat_output_dir):
        taxa_dir_path = os.path.join(master_agat_output_dir, taxa_dir)

        for gene_dir in listdir(taxa_dir_path):
            gene_path = os.path.join(taxa_dir_path, gene_dir)

            for file in listdir(gene_path):

                filepath = os.path.join(gene_path, file)
                f = open(filepath, "r")

                line = f.readline()
                if len(line) == 0 or line[0] != ">":
                    print("Empty or ill-formatted file: " + gene_dir + " " + file)
                    print("Continuing as normal...")

                else:

                    combined_filepath = os.path.join(combined_agat_file_dir, gene_dir + ".fas")
                    f2 = open(combined_filepath, "a")

                    while line != "":

                        fasta_heading = ">" + os.path.splitext(file)[0] + " " + line.strip()[1:]
                        sequence = ""

                        line = f.readline()
                        while line != "" and line[0] != ">":
                            if line.strip() != "":
                                sequence += line.strip()
                            line = f.readline()

                        f2.write(fasta_heading + "\n" + sequence + "\n")

    return combined_agat_file_dir


def add_references_to_combined_files(reference_sequences_fasta, combined_agat_file_dir):

    """
    Add reference sequences to .fasta files containing gene predictions grouped by gene.

    :param reference_sequences_fasta: Directory containing reference sequences
    :param combined_agat_file_dir: Directory containing combined predictions by gene.
    :return: Output directory containing combined predictions by gene, with reference sequences added.
    """

    for reference_taxon_dir in listdir(reference_sequences_fasta):
        reference_taxon_path = os.path.join(reference_sequences_fasta, reference_taxon_dir)

        for reference_gene_file in listdir(reference_taxon_path):
            reference_gene_path = os.path.join(reference_taxon_path, reference_gene_file)

            if reference_gene_file in listdir(combined_agat_file_dir):
                existing_combined_filestream = open(os.path.join(combined_agat_file_dir, reference_gene_file), "a")

                reference_tagged_stream = ">REFERENCE " + open(reference_gene_path, "r").read()[1:]
                existing_combined_filestream.write(reference_tagged_stream)


def seq_filter(mode,
               reference_sequences_fasta,
               combined_agat_file_dir,
               download_dir,
               score_threshold=1.0,
               lower_length_bound=1, upper_length_bound=1):
    """
    Filter out poor-quality gene predictions. Here, quality is defined by length similarity (comparing the lengths of
    the prediction to the reference), or by sequence similarity (globally aligning the prediction against the reference
    and observing the resultant score).

    :param mode: Either "length" (number of bases) or "score" (raw score from global alignment)
    :param score_threshold: Filter threshold used if "score" is selected as the mode. All gene predictions with scores
    below the score_threshold multiplied by the maximum possible raw score (when the reference sequence is used) are
    filtered out.
    :param upper_length_bound: Upper filter threshold if "length" is selected as the mode. All gene predictions with
    lengths above the reference sequence length multiplied by upper_length_bound are removed.
    :param lower_length_bound: Lower filter threshold if "length" is selected as the mode. All gene predictions with
    lengths below the reference sequence length multiplied by upper_length_bound are removed.
    :param reference_sequences_fasta: Directory containing reference sequences in nucleotide format. The master
    directory must contain separate child directories for each taxa, and each child directory must contain only .fasta
    files.
    :param combined_agat_file_dir: Directory containing combined gene predictions as nucleotide sequences.
    :param download_dir: Path to download output.
    :return: Directory containing gene prediction fasta files with poor quality sequences removed.
    """
    # create a new directory

    if mode == "length":
        master_filtered_agat_dir = (
            make_unique_directory(os.path.join(download_dir, "simple_length_filter_" + str(lower_length_bound) + "_" + str(upper_length_bound))))
    else:
        master_filtered_agat_dir = (
            make_unique_directory(os.path.join(download_dir, "score_filter_" + str(score_threshold))))

    for combined_agat_file in listdir(combined_agat_file_dir):
        combined_agat_file_path = os.path.join(combined_agat_file_dir, combined_agat_file)

        reference_seqs= []

        for reference_taxon_dir in listdir(reference_sequences_fasta):
            reference_taxon_path = os.path.join(reference_sequences_fasta, reference_taxon_dir)

            reference_gene_files = listdir(reference_taxon_path)
            if combined_agat_file in reference_gene_files:

                reference_stream = open(os.path.join(reference_taxon_path, combined_agat_file), "r")
                reference_seq = reference_stream.readlines()[1].strip()
                reference_seqs.append(reference_seq)

        filtered_agat_file = os.path.join(master_filtered_agat_dir, combined_agat_file)
        filtered_agat_stream = open(filtered_agat_file, "a")

        combined_file_arr = open(combined_agat_file_path, "r").readlines()
        for i in range(1, len(combined_file_arr), 2):

            stripped_seq = combined_file_arr[i].strip()

            passed = True
            for r_seq in reference_seqs:

                if mode == "length":
                    if len(stripped_seq) < len(r_seq) * lower_length_bound or  len(stripped_seq) > len(r_seq) * upper_length_bound:
                        passed = False
                        break

                else:
                    aligner = Bio.Align.PairwiseAligner()
                    aligner.match_score = 2
                    aligner.mismatch_score = -3
                    aligner.open_gap_score = -5
                    aligner.extend_gap_score = -2

                    score = aligner.score(r_seq, stripped_seq)
                    if score < aligner.match_score * len(r_seq) * score_threshold:
                        passed = False
                        break

            if passed:
                filtered_agat_stream.write(combined_file_arr[i-1] + combined_file_arr[i])

    return master_filtered_agat_dir



if __name__ == "__main__":
    """
    seq_filter(1, 1,
                         "/Users/tonyx/Documents/chang_lab/miniprot_reference_nuc_with_variants",
                         "/Users/tonyx/Documents/chang_lab/combinated_agat_files_variants",
                         "/Users/tonyx/Documents/chang_lab")
    """

    seq_filter("score", "/Users/tonyx/Documents/chang_lab/miniprot_reference_nuc_with_variants",
                         "/Users/tonyx/Documents/chang_lab/combinated_agat_files_variants",
                         "/Users/tonyx/Documents/chang_lab", score_threshold=0.4)












