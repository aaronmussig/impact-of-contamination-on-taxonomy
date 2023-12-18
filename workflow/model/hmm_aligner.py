
from Bio import AlignIO
import io
def get_aligned_marker(hit_name, result_file):
    """
    Parse the output of Hmmalign
    :param hit_name: gene name
    :param result_file: output file from Hmmalign
    """
    hit_seq = None
    mask_seq = None

    for line in result_file.splitlines():
        splitline = line.split(" ", 1)
        if splitline[0] == hit_name.split(" ", 1)[0]:
            rsplitline = line.rsplit(" ", 1)
            hit_seq = rsplitline[-1]
        elif line[0:len("#=GC RF")] == "#=GC RF":
            rsplitline = line.rsplit(" ", 1)
            mask_seq = rsplitline[-1]

    if mask_seq is None:
        raise Exception("Unable to get mask from hmm align result file")

    if hit_seq is None:
        return None

    aligned_marker = ''.join([h for h, m in zip(hit_seq, mask_seq) if m == 'x'])
    return aligned_marker

def get_aligned_markers(result_file):
    """
    Parse the output of Hmmalign
    :param hit_name: gene name
    :param result_file: output file from Hmmalign
    """
    result_io = io.StringIO(result_file)
    align = AlignIO.read(result_io, "stockholm")

    out = dict()
    for record in align:
        seq_out = list()
        for aa, mask in zip(str(record.seq), align.column_annotations['reference_annotation']):
            if mask == 'x':
                seq_out.append(aa)
        out[record.id] = ''.join(seq_out)
    return out
