import os
import sys
from typing import Optional
#
# # --------------------------------------------------------------------------------
# # Credentials
# # --------------------------------------------------------------------------------
#
REDIS_HOST = os.environ.get('REDIS_HOST')
REDIS_PASS = os.environ.get('REDIS_PASS')
#
# # --------------------------------------------------------------------------------
# # Environment variables
# # --------------------------------------------------------------------------------
#
# # The remote env name of the gunc chimeras
REMOTE_ENV = 'gunc-chimeras'
#
# CONDA_PATH = '/opt/centos7/sw/miniconda3/etc/profile.d/conda.sh'
# if not os.path.exists(CONDA_PATH):
#     CONDA_PATH = '/Users/aaron/opt/miniconda3/etc/profile.d/conda.sh'
#
# --------------------------------------------------------------------------------
# Project path variables
# --------------------------------------------------------------------------------

# Define the project root directory
DIR_PROJECT = '/srv/home/uqamussi/projects/gunc-chimeras'

DIR_OUT_ROOT = '/srv/home/uqamussi/projects/gunc-chimeras/output'
DIR_OUT_EXTERNAL = os.path.join(DIR_OUT_ROOT, 'external')
DIR_OUT_GTDB = os.path.join(DIR_OUT_ROOT, 'gtdb')
DIR_OUT_BATCH = os.path.join(DIR_OUT_ROOT, 'batch')
DIR_OUT_GUNC = os.path.join(DIR_OUT_ROOT, 'gunc')
DIR_OUT_CUTOFF = os.path.join(DIR_OUT_ROOT, 'cutoff')
DIR_OUT_FASTANI = os.path.join(DIR_OUT_ROOT, 'fastani')
DIR_OUT_FASTANI_INTER = os.path.join(DIR_OUT_ROOT, 'fastani_interspecies')
DIR_OUT_HMMER = os.path.join(DIR_OUT_ROOT, 'hmmer')
DIR_OUT_MARKER = os.path.join(DIR_OUT_ROOT, 'marker')
DIR_OUT_TREE = os.path.join(DIR_OUT_ROOT, 'tree')
DIR_OUT_PROGENOMES = os.path.join(DIR_OUT_ROOT, 'progenomes')
DIR_OUT_BACKGROUND = os.path.join(DIR_OUT_ROOT, 'background')
DIR_OUT_BOOTSTRAP = os.path.join(DIR_OUT_ROOT, 'bootstrap')

# The directory where the source code is located (accessible by RQ workers).
DIR_CODE = os.path.join(DIR_PROJECT, 'code')

DIR_OUT_SENTINEL = os.path.join(DIR_OUT_ROOT, 'sentinel')


PCT_VALUES = (1, 5, 10, 15, 20, 30, 40, 50)

# # Directory to store associated luigi files
# DIR_LUIGI = os.path.join(DIR_ROOT, 'luigi')
#
# # Directory where the root genomes and working data will be stored
DIR_OUT_GENOMES = os.path.join(DIR_OUT_ROOT, 'genomes')
#
# # This directory is used for external files (i.e. files that can be downloaded again)
# DIR_EXTERNAL = os.path.join(DIR_ROOT, 'external')
#
# # The directory that will contain all NCBI genomes
# DIR_EXTERNAL_GENOMES = os.path.join(DIR_EXTERNAL, 'genomes')
#
# # Only used to create output files for jobs that impact large number of files
# DIR_LUIGI_SENTINEL = os.path.join(DIR_LUIGI, 'sentinel')
# DIR_LUIGI_OUTPUT = os.path.join(DIR_LUIGI, 'output')
# DIR_LUIGI_ANALYSIS = os.path.join(DIR_LUIGI, 'analysis')
# DIR_LUIGI_INTERMEDIATE = os.path.join(DIR_LUIGI, 'intermediate')
# DIR_LUIGI_PRODIGAL = os.path.join(DIR_LUIGI, 'prodigal')
#
# # If this is set, then files will be copied locally to this directory on access
DIR_CACHE: Optional[str] = '/tmp/gunc-cache'
# DIR_CACHE = None

#
# True if the code is being run through a debugger
DEBUG = sys.gettrace() is not None
if DEBUG:
    os.nice(19)
#
# # --------------------------------------------------------------------------------
# # Workflow variables
# # --------------------------------------------------------------------------------
#
# GUNC_BATCH_SIZE = 20000
# DIR_GUNC_BATCH = os.path.join(DIR_LUIGI_INTERMEDIATE, 'gunc_batches')
# GUNC_CUTOFF_BATCH_DIR = os.path.join(DIR_LUIGI_INTERMEDIATE, 'gunc_cutoff_batches')
# DIR_GUNC_CUTOFF_H5 = os.path.join(DIR_LUIGI_OUTPUT, 'gunc_cutoff')
#
# # --------------------------------------------------------------------------------
# # ACE-specific variables
# # --------------------------------------------------------------------------------
#
# True if you've got access to the ACE file share.
AT_ACE = True
#
# # Path to the ACE file share NCBI mirror
ACE_R207_ROOT = '/srv/db/gtdb/genomes/ncbi/release207'
ACE_R95_ROOT = '/srv/db/gtdb/genomes/ncbi/release95'
NCBI_FTP_ROOT = '/srv/db/ncbi/new_ftp_structure/genomes/all'

#
# # --------------------------------------------------------------------------------
# # Reference databases
# # --------------------------------------------------------------------------------
#
# # Path to the GUNC diamond reference database.
PATH_GUNC_GTDB_REF_DB = '/srv/db/gunc/gunc_db_gtdb95.dmnd'
PATH_GUNC_PRO_REF_DB = '/srv/db/gunc/gunc_db_progenomes2.1.dmnd'

#
# # --------------------------------------------------------------------------------
# # RQ-specific variables
# # --------------------------------------------------------------------------------
#
# # Time to wait until checking if the Job is still running in RQ
RQ_SLEEP_DELAY = 120
#
# RQ_GENERAL = 'general'
#
# RQ_CPU_1 = 'cpu-1'  # single threaded, will auto-scale up on each server
# RQ_CPU_N = 'cpu-n'  # multithreaded, will only execute one job per server
# RQ_MEM_HIGH = 'mem-high'  # high memory jobs, will only execute one job per server if memory is available
#
# # Path to GUNC run on GTDB R95
# DIR_GUNC_R95 = os.path.join(DIR_ROOT, 'gtdb_r95')
#
# DIR_GUNC_R95_OLD = '/srv/home/uqamussi/projects/gunc-chimeras/gtdb_r95'
#
# # Redis Queue
# RQ_FASTANI_AFTER_REMOVING_CONTIGS = 'fastani_after_removing_contigs'
# RQ_HMM_AFTER_REMOVING_CONTIGS = 'gunc_hmm_remove_contigs_2'
#
R95_AF_THRESHOLD = 0.65
R207_AF_THRESHOLD = 0.5

#
# # Contains additional genomes that were missing in the release95 directory
# DIR_R95_FIX = os.path.join(DIR_ROOT, 'release95')
#
R95_ROOT = '/srv/db/gtdb/genomes/ncbi/release95'
R202_ROOT = '/srv/db/gtdb/genomes/ncbi/release207'
#
# # External static files
# DIR_EXTERNAL = os.path.join(DIR_ROOT, 'external_data')
#
# DIR_EXTERNAL_GTDB = os.path.join(DIR_EXTERNAL, 'gtdb')
#
# # Paths to databases for Pfam and Tigrfam
DIR_PFAM_27 = '/srv/db/pfam/27'
DIR_PFAM_33 = '/srv/db/pfam/33.1'
DIR_TIGRFAM_15 = '/srv/db/tigrfam/15.0/TIGRFAMs_15.0_HMM'

# # Set in the contig assignment file, used as indices
TAX_LEVELS = {'kingdom': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6}

R95_BAC120_HMM = {'PF00380.14': 121, 'PF00410.14': 129, 'PF00466.15': 100, 'PF01025.14': 166, 'PF02576.12': 141,
                  'PF03726.9': 83, 'TIGR00006': 310, 'TIGR00019': 361, 'TIGR00020': 365, 'TIGR00029': 87,
                  'TIGR00043': 111, 'TIGR00054': 421, 'TIGR00059': 112, 'TIGR00061': 101, 'TIGR00064': 279,
                  'TIGR00065': 353, 'TIGR00082': 115, 'TIGR00083': 290, 'TIGR00084': 192, 'TIGR00086': 144,
                  'TIGR00088': 233, 'TIGR00090': 99, 'TIGR00092': 368, 'TIGR00095': 194, 'TIGR00115': 410,
                  'TIGR00116': 293, 'TIGR00138': 183, 'TIGR00158': 148, 'TIGR00166': 95, 'TIGR00168': 165,
                  'TIGR00186': 240, 'TIGR00194': 574, 'TIGR00250': 130, 'TIGR00337': 526, 'TIGR00344': 847,
                  'TIGR00362': 437, 'TIGR00382': 414, 'TIGR00392': 861, 'TIGR00396': 843, 'TIGR00398': 530,
                  'TIGR00414': 418, 'TIGR00416': 454, 'TIGR00420': 351, 'TIGR00431': 210, 'TIGR00435': 466,
                  'TIGR00436': 270, 'TIGR00442': 406, 'TIGR00445': 321, 'TIGR00456': 569, 'TIGR00459': 586,
                  'TIGR00460': 315, 'TIGR00468': 324, 'TIGR00472': 798, 'TIGR00487': 587, 'TIGR00496': 176,
                  'TIGR00539': 361, 'TIGR00580': 923, 'TIGR00593': 890, 'TIGR00615': 196, 'TIGR00631': 658,
                  'TIGR00634': 563, 'TIGR00635': 305, 'TIGR00643': 629, 'TIGR00663': 367, 'TIGR00717': 516,
                  'TIGR00755': 256, 'TIGR00810': 73, 'TIGR00922': 172, 'TIGR00928': 436, 'TIGR00959': 428,
                  'TIGR00963': 787, 'TIGR00964': 57, 'TIGR00967': 414, 'TIGR01009': 212, 'TIGR01011': 225,
                  'TIGR01017': 200, 'TIGR01021': 156, 'TIGR01029': 154, 'TIGR01032': 114, 'TIGR01039': 462,
                  'TIGR01044': 103, 'TIGR01059': 639, 'TIGR01063': 800, 'TIGR01066': 141, 'TIGR01071': 144,
                  'TIGR01079': 104, 'TIGR01082': 449, 'TIGR01087': 441, 'TIGR01128': 314, 'TIGR01146': 286,
                  'TIGR01164': 126, 'TIGR01169': 227, 'TIGR01171': 275, 'TIGR01302': 450, 'TIGR01391': 414,
                  'TIGR01393': 595, 'TIGR01394': 594, 'TIGR01510': 155, 'TIGR01632': 140, 'TIGR01951': 131,
                  'TIGR01953': 340, 'TIGR02012': 321, 'TIGR02013': 1238, 'TIGR02027': 298, 'TIGR02075': 233,
                  'TIGR02191': 219, 'TIGR02273': 166, 'TIGR02350': 596, 'TIGR02386': 1147, 'TIGR02397': 355,
                  'TIGR02432': 189, 'TIGR02729': 329, 'TIGR03263': 180, 'TIGR03594': 432, 'TIGR03625': 202,
                  'TIGR03632': 117, 'TIGR03654': 175, 'TIGR03723': 314, 'TIGR03725': 212, 'TIGR03953': 188}
R95_AR122_HMM = {'PF00368.13': 373, 'PF00410.14': 129, 'PF00466.15': 100, 'PF00687.16': 220, 'PF00827.12': 192,
                 'PF00900.15': 77, 'PF01000.21': 112, 'PF01015.13': 195, 'PF01090.14': 140, 'PF01092.14': 127,
                 'PF01157.13': 99, 'PF01191.14': 74, 'PF01194.12': 60, 'PF01198.14': 83, 'PF01200.13': 69,
                 'PF01269.12': 229, 'PF01280.15': 148, 'PF01282.14': 84, 'PF01496.14': 759, 'PF01655.13': 110,
                 'PF01798.13': 150, 'PF01864.12': 175, 'PF01866.12': 307, 'PF01868.11': 89, 'PF01984.15': 107,
                 'PF01990.12': 95, 'PF02006.11': 178, 'PF02978.14': 104, 'PF03874.11': 117, 'PF04019.7': 121,
                 'PF04104.9': 260, 'PF04919.7': 181, 'PF07541.7': 114, 'PF13656.1': 77, 'PF13685.1': 250,
                 'TIGR00021': 218, 'TIGR00037': 130, 'TIGR00042': 184, 'TIGR00064': 279, 'TIGR00111': 351,
                 'TIGR00134': 622, 'TIGR00240': 150, 'TIGR00264': 116, 'TIGR00270': 154, 'TIGR00279': 172,
                 'TIGR00283': 115, 'TIGR00291': 231, 'TIGR00293': 129, 'TIGR00307': 127, 'TIGR00308': 375,
                 'TIGR00323': 215, 'TIGR00324': 177, 'TIGR00335': 324, 'TIGR00336': 173, 'TIGR00337': 526,
                 'TIGR00373': 162, 'TIGR00389': 565, 'TIGR00392': 861, 'TIGR00398': 530, 'TIGR00405': 145,
                 'TIGR00408': 475, 'TIGR00422': 863, 'TIGR00425': 322, 'TIGR00432': 637, 'TIGR00442': 406,
                 'TIGR00448': 179, 'TIGR00456': 569, 'TIGR00458': 428, 'TIGR00463': 560, 'TIGR00468': 324,
                 'TIGR00471': 551, 'TIGR00490': 720, 'TIGR00491': 594, 'TIGR00501': 295, 'TIGR00521': 392,
                 'TIGR00522': 258, 'TIGR00549': 276, 'TIGR00658': 304, 'TIGR00670': 304, 'TIGR00729': 207,
                 'TIGR00936': 416, 'TIGR00982': 139, 'TIGR01008': 195, 'TIGR01012': 196, 'TIGR01018': 162,
                 'TIGR01020': 212, 'TIGR01025': 135, 'TIGR01028': 186, 'TIGR01038': 148, 'TIGR01046': 99,
                 'TIGR01052': 488, 'TIGR01060': 425, 'TIGR01077': 141, 'TIGR01080': 116, 'TIGR01213': 387,
                 'TIGR01309': 151, 'TIGR01952': 141, 'TIGR02076': 222, 'TIGR02153': 405, 'TIGR02236': 311,
                 'TIGR02258': 180, 'TIGR02338': 110, 'TIGR02389': 367, 'TIGR02390': 867, 'TIGR02651': 302,
                 'TIGR03626': 331, 'TIGR03627': 130, 'TIGR03628': 117, 'TIGR03629': 144, 'TIGR03636': 77,
                 'TIGR03653': 170, 'TIGR03665': 173, 'TIGR03670': 599, 'TIGR03671': 410, 'TIGR03672': 251,
                 'TIGR03673': 131, 'TIGR03674': 338, 'TIGR03677': 117, 'TIGR03680': 407, 'TIGR03683': 902,
                 'TIGR03684': 154, 'TIGR03722': 323}
assert (len(R95_AR122_HMM) == 122)
assert (len(R95_BAC120_HMM) == 120)
