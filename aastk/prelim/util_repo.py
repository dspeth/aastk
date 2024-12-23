def get_num_sequences(fasta_file):
    """Return the number of sequences in a fasta file."""
    return sum(1 for line in open(fasta_file) if line.startswith(">"))

def check_tool_installed(tool_name):
    """Check whether the tool (e.g., 'diamond', 'blastp') is available in the PATH."""
    if which(tool_name) is None:
        raise EnvironmentError(f"{tool_name} is not available in the system PATH")

def validate_file(file_path, file_type="fasta"):
    """Validate if a file exists and is of the correct type."""
    file = Path(file_path)
    if not file.is_file():
        raise FileNotFoundError(f"{file_path} does not exist")
    if file_type == "fasta":
        with open(file) as f:
            if not f.readline().startswith(">"):
                raise ValueError(f"{file_path} does not seem to be a FASTA file")


def fasta_subsample(fasta_file, output_file, subset_size):
    """Randomly sample 'subset_size' sequences from a fasta file and write to output_file."""
    seq_index = SeqIO.index(fasta_file, "fasta")
    subset_keys = random.sample(seq_index.keys(), subset_size)
    subset_dict = {k: seq_index[k] for k in subset_keys}
    SeqIO.write(subset_dict.values(), output_file, "fasta-2line")