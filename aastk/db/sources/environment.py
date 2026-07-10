import logging
logger = logging.getLogger(__name__)

def parse_environment_file(env_filepath: str, columns: list) -> list:
    environment_data = []
    with open(env_filepath, 'r') as f:
        next(f)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(columns) + 1:
                continue

            genome_id = parts[0]
            row = [genome_id]
            for i, col in enumerate(columns, start=1):
                value = float(parts[i]) if parts[i].strip() else None
                row.append(value)

            environment_data.append(tuple(row))

    return environment_data
