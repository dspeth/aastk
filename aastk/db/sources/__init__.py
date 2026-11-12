QUOTE_TRANSLATION = str.maketrans({'‘': "'", '’': "'", '“': '"', '”': '"'})


def read_lines_utf8_or_cp1252(filepath: str) -> list:
    with open(filepath, 'rb') as f:
        raw = f.read()
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        text = raw.decode('cp1252', errors='replace')
    return text.splitlines()
