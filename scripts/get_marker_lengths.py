from workflow.config import R207_MARKERS, DIR_PFAM_33, DIR_TIGRFAM_15
import os

def _get_hmm_size(path):
    size = 0
    with open(path) as fp:
        for line in fp:
            if line.startswith("LENG  "):
                size = line.split("  ")[1]
                break
    return int(size)



def main():

    print("R207_MARKER_LENGTHS = {")

    for marker in R207_MARKERS:
        if marker.startswith('PF'):
            hmm_path = os.path.join(DIR_PFAM_33, 'individual_hmms', f'{marker}.hmm')
        elif marker.startswith('TIGR'):
            hmm_path = os.path.join(DIR_TIGRFAM_15, 'individual_hmms', f'{marker}.HMM')
        else:
            raise Exception('Unknown marker')

        size = _get_hmm_size(hmm_path)
        print(f"    '{marker}': {size},")

    print("}")
    return

if __name__ == '__main__':
    main()