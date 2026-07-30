# input: protein fasta
# optional arguments: output dir, model, batch size (sequences) depending on sequence length?, force/keep
# output: .h5 (real embedding output), metadata json?, tsv for quick inspection?

from pathlib import Path

import torch
import h5py
from transformers import AutoTokenizer, T5EncoderModel

from aastk.util import ensure_path, stream_fasta, clean_fasta_header

import logging
logger = logging.getLogger(__name__)

MODEL_REGISTRY = {
    "prott5": "Rostlab/prot_t5_xl_half_uniref50-enc"
}

SEQUENCE_LENGTH_WARNING = 3000

def resolve_output_path(
        fasta_path: Path,
        output: str = None,
        force: bool = False):

    fasta_name = fasta_path.name

    if fasta_name.endswith(".gz"):
        fasta_name = Path(fasta_name).stem

    fasta_stem = Path(fasta_name).stem
    output_file = f"{fasta_stem}_embeddings.h5"

    # no output provided as argument
    if output is None:
        return Path(ensure_path(
            path=str(fasta_path.parent),
            target=output_file,
            force=force
        ))

    output_path = Path(output)

    # output file name is provided
    if output_path.suffix.lower() in [".h5", ".hdf5"]:
        return Path(ensure_path(
            path=str(output_path.parent),
            target=output_path.name,
            force=force
        ))

    # output is something else -> handled as directory
    return Path(ensure_path(
            path=str(output_path),
            target=output_file,
            force=force
    ))

def prepare_sequence_prott5(sequence: str):
    """
    prott5 expects amino acids separated by spaces
    amino acids U, Z and O are mapped to X
    """

    sequence = "".join(sequence.split()).upper().replace("-", "")
    sequence = sequence.replace("U", "X").replace("Z", "X").replace("O", "X")

    return sequence

def add_spaces(sequence: str):
    """
    adds spaces between residues to fit ProtT5 formatting
    """
    return " ".join(sequence)

def h5_safe_id(protein_id: str):
    """
    HDF5 object names must not contain slashes
    """

    return protein_id.replace("/", "_slash_")

def collect_fasta_records(fasta_path: Path):
    """
    reading protein FASTA records and preparing them
    returned list-entries are a 3-tuple: protein_id, tokenized_sequence, sequence_length
    """

    if not fasta_path.is_file():
        raise FileNotFoundError(f"Input FASTA does not exist: {fasta_path}")

    records = []
    seen_ids = set()

    for header, sequence in stream_fasta(str(fasta_path)):
        protein_id = h5_safe_id(clean_fasta_header(header))

        cleaned_sequence = prepare_sequence_prott5(sequence)
        sequence_length = len(cleaned_sequence)
        model_ready_sequence = add_spaces(cleaned_sequence)

        if not protein_id:
            raise ValueError("FASTA contains an empty header")

        if protein_id in seen_ids:
            raise ValueError(f"Duplicate protein ID in FASTA: {protein_id}")

        if not cleaned_sequence:
            raise ValueError(f"FASTA contains an empty sequence: {protein_id}")

        if sequence_length > SEQUENCE_LENGTH_WARNING:
            logger.warning(
                f"Sequence of {protein_id} is longer than the maximum sequence length: {SEQUENCE_LENGTH_WARNING}"
            )

        seen_ids.add(protein_id)
        records.append((protein_id, model_ready_sequence, sequence_length))

    if not records:
        raise ValueError(f"FASTA contains no sequences: {fasta_path}")

    logger.info(f"Read {len(records)} protein sequence(s) from {fasta_path}")

    return records


def load_plm_model(model: str):
    """
    loading PLM model and respective Tokenizer
    """

    model_name = MODEL_REGISTRY[model]

    logger.info(f"Loading PLM model: {model_name}")

    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    plm_model = T5EncoderModel.from_pretrained(model_name)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if device.type == "cpu":
        plm_model = plm_model.float()

    plm_model = plm_model.to(device)
    plm_model.eval()

    logger.info(f"Using device: {device}")

    return plm_model, tokenizer, device

def batch_handling(records: list[tuple[str, str, int]], max_residues_per_batch: int):
    """
    creating batches while limiting sequence/residue count
    """

    if max_residues_per_batch < 1:
        raise ValueError("max_residues_per_batch must be at least 1 (default: 3000)")

    current_batch = []
    current_residues = 0

    # sort, long sequences first
    records = sorted(records, key=lambda record: record[2], reverse=True)

    for record in records:
        sequence_length = record[2]

        # if sequence is longer than residue limit -> process alone
        if sequence_length > max_residues_per_batch:
            if current_batch:
                yield current_batch
                current_batch = []
                current_residues = 0

            logger.warning(f"Sequence {record[0]} has {sequence_length} residues and exceeds the limit of {max_residues_per_batch}. Processing it alone")

            yield [record]
            continue

        residue_limit_reached = (
                current_batch
                and current_residues + sequence_length > max_residues_per_batch
        )

        if residue_limit_reached:
            yield current_batch
            current_batch = []
            current_residues = 0

        current_batch.append(record)
        current_residues += sequence_length

    if current_batch:
        yield current_batch

def embed_current_batch(batch, plm_model, tokenizer, device):
    """
    do model call for current batch and return embeddings
    """

    protein_ids, sequences, sequence_lengths = zip(*batch)

    token_encoding = tokenizer(
        list(sequences),
        add_special_tokens=True,
        padding="longest",
        return_tensors="pt",
    )

    input_ids = token_encoding["input_ids"].to(device)
    attention_mask = token_encoding["attention_mask"].to(device)

    with torch.no_grad():
        model_output = plm_model(input_ids=input_ids, attention_mask=attention_mask)

    embeddings = {}

    for batch_index, protein_id in enumerate(protein_ids):
        sequence_length = sequence_lengths[batch_index]

        # Slice away padding and the special token, then mean-pool over residues.
        residue_embeddings = model_output.last_hidden_state[batch_index, :sequence_length]
        mean_embedding = residue_embeddings.mean(dim=0)

        embeddings[protein_id] = mean_embedding.detach().cpu().numpy().astype("float32")

    return embeddings

def write_embeddings(records, output_path: Path, plm_model, tokenizer, device, max_residues_per_batch: int):
    """
    writing mean-pooled protein embeddings to h5
    repeatedly call embed_current_batch
    """

    embedding_count = 0

    with h5py.File(output_path, "w") as h5_file:
        batches = batch_handling(records, max_residues_per_batch)

        for batch_number, batch in enumerate(batches, start=1):
            batch_residues = sum(record[2] for record in batch)

            logger.info(f"Embedding batch {batch_number} with {len(batch)} sequence(s) and {batch_residues} residue(s)")

            embeddings = embed_current_batch(
                batch=batch,
                plm_model=plm_model,
                tokenizer=tokenizer,
                device=device,
            )

            for protein_id, embedding in embeddings.items():
                h5_file.create_dataset(protein_id, data=embedding)
                embedding_count += 1

    logger.info(f"Wrote {embedding_count} embedding(s) to {output_path}")
    return output_path

def plm_embedder(
        fasta: str,
        model: str,
        output: str = None,
        max_residues_per_batch: int = 3000,
        force: bool = False):
    logger.info("Starting PLM embedder")

    fasta_path = Path(fasta)
    output_path = resolve_output_path(fasta_path=fasta_path, output = output, force=force)

    records = collect_fasta_records(fasta_path=fasta_path)

    plm_model, tokenizer, device = load_plm_model(model=model)

    write_embeddings(
        records=records,
        output_path=output_path,
        plm_model=plm_model,
        tokenizer=tokenizer,
        device=device,
        max_residues_per_batch=max_residues_per_batch
    )

    logger.info("PLM embedder complete")

    return str(output_path)

