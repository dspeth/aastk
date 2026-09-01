# define database table schema
PROTEIN_SCHEMA = {
    'seqID': 'TEXT PRIMARY KEY',
    'parent_ID': 'TEXT',
    'aa_length': 'INTEGER',
    'strand': 'TEXT',
    'COG_ID': 'TEXT',
    'KEGG_ID': 'TEXT',
    'Pfam_ID': 'TEXT',
    'cugo_number': 'INTEGER',
    'no_tmh': 'INTEGER',
    'protein_seq': 'BLOB'
}

GENOME_SCHEMA = {
    'genome_ID': 'TEXT PRIMARY KEY',
    'domain': 'TEXT',
    'phylum': 'TEXT',
    'class': 'TEXT',
    'order_tax': 'TEXT',
    'family': 'TEXT',
    'genus': 'TEXT',
    'species': 'TEXT',
    'culture_collection': 'INT',
    'animal_associated': 'REAL',
    'aquatic': 'REAL',
    'built': 'REAL',
    'other': 'REAL',
    'sediment': 'REAL',
    'soil': 'REAL',
    'unassigned_high_level': 'REAL',
    'human': 'REAL',
    'invertebrate': 'REAL',
    'other_vertebrate': 'REAL',
    'unspecified_animal': 'REAL',
    'aquatic_other': 'REAL',
    'aquatic_unspecified': 'REAL',
    'freshwater': 'REAL',
    'groundwater': 'REAL',
    'marine': 'REAL',
    'built_other': 'REAL',
    'drinking_water': 'REAL',
    'wastewater': 'REAL',
    'air': 'REAL',
    'bacteria': 'REAL',
    'eukaryote_other': 'REAL',
    'food': 'REAL',
    'geothermal': 'REAL',
    'hypersaline': 'REAL',
    'other_unspecified': 'REAL',
    'plant_associated': 'REAL',
    'subsurface': 'REAL',
    'synthetic': 'REAL',
    'viral': 'REAL',
    'freshwater_sediment': 'REAL',
    'marine_sediment': 'REAL',
    'sediment_unspecified': 'REAL',
    'desert': 'REAL',
    'forest': 'REAL',
    'rhizosphere': 'REAL',
    'soil_agricultural': 'REAL',
    'soil_other': 'REAL',
    'soil_unspecified': 'REAL',
    'tundra_wetland': 'REAL',
    'unassigned_low_level': 'REAL'
}

# extract column name lists for queries (same order)
PROTEIN_COLUMNS = list(PROTEIN_SCHEMA.keys())
GENOME_COLUMNS = list(GENOME_SCHEMA.keys())

# categorize columns for later use
BASE_COLUMNS = ['seqID', 'parent_ID', 'aa_length', 'strand']

ANNOTATION_COLUMNS = ['COG_ID', 'KEGG_ID', 'Pfam_ID']

TAXONOMY_COLUMNS = [
    'genome_ID', 'domain', 'phylum', 'class',
    'order_tax', 'family', 'genus', 'species'
]

CULTURE_COLLECTION_COLUMNS = ['culture_collection']

HIGH_LEVEL_ENV_COLUMNS = [
    'animal_associated', 'aquatic', 'built', 'other',
    'sediment', 'soil', 'unassigned_high_level'
]

LOW_LEVEL_ENV_COLUMNS = [
    'human', 'invertebrate', 'other_vertebrate', 'unspecified_animal',
    'aquatic_other', 'aquatic_unspecified', 'freshwater', 'groundwater', 'marine',
    'built_other', 'drinking_water', 'wastewater',
    'air', 'bacteria', 'eukaryote_other', 'food', 'geothermal', 'hypersaline',
    'other_unspecified', 'plant_associated', 'subsurface', 'synthetic', 'viral',
    'freshwater_sediment', 'marine_sediment', 'sediment_unspecified',
    'desert', 'forest', 'rhizosphere', 'soil_agricultural', 'soil_other',
    'soil_unspecified', 'tundra_wetland', 'unassigned_low_level'
]

HIGH_LEVEL_ENV_VIEW = 'genome_high_level_env'
LOW_LEVEL_ENV_VIEW = 'genome_low_level_env'
HIGH_LEVEL_ENV_CATEGORY_COLUMN = 'high_level_env_category'
LOW_LEVEL_ENV_CATEGORY_COLUMN = 'low_level_env_category'
