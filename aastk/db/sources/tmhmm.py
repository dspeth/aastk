import logging

logger = logging.getLogger(__name__)

def process_tmhmm_file(tmhmm_filepath):
    tmhmm_data = []
    try:
        with open(tmhmm_filepath, 'r') as f:
            for line_no, line in enumerate(f, start=1):
                if line.startswith('prot_ID'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    prot_id = parts[0]
                    try:
                        no_tmh = int(parts[1])
                    except ValueError:
                        logger.warning(
                            f"{tmhmm_filepath}:{line_no}: non-integer no_tmh value {parts[1]!r}, skipping line"
                        )
                        continue
                    tmhmm_data.append((no_tmh, prot_id))
    except OSError as e:
        logger.warning(f"Failed to read TMHMM file {tmhmm_filepath}: {e}")
        return []
    return tmhmm_data