#!/usr/bin/env python3

if __name__ == '__main__':
    import sys

    if len(sys.argv) != 5:
        sys.exit('Usage: coverage2bed.py <infile> <depth> <minfrac> <outfile>')

    minfrac = float(sys.argv[3])

    state = 0
    current_chrom = None
    last_pos = None

    with open(sys.argv[1], 'rt') as infile:
        header = next(infile).rstrip().split('\t')
        assert sys.argv[2] in header, 'Depth "{}" not found in coverage file header'.format(sys.argv[2])
        frac_idx = header.index(sys.argv[2])
        with open(sys.argv[4], 'wt') as outfile:
            for line in infile:
                fields = line.rstrip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                frac = float(fields[frac_idx])

                # State machine
                if chrom != current_chrom:
                    # Changed chromosomes.  Transition to state 0, to be picked up by the 
                    # following code in this same loop iteration.
                    if state == 1:
                        outfile.write('{}\n'.format(last_pos))
                    state = 0
                    current_chrom = chrom

                if state == 0:
                    # Initial state at the start of a chromosome. Decide whether 
                    # to transition to state 1 or 2, depending on frac.
                    if frac >= minfrac:
                        # Move to state 1 (start of an output bed region). Emit the start of the bed line.
                        outfile.write('{}\t{}\t'.format(chrom, pos - 1))
                        state = 1
                    else:
                        # Move to state 2 (not in an output bed region).
                        state = 2
                elif state == 1:
                    # Inside a valid region. Verify that this is still the case.
                    if frac < minfrac:
                        # Outside of a valid region now, but still on the same chrom.
                        # Emit the end of the region and transition to state 2.
                        outfile.write('{}\n'.format(last_pos))
                        state = 2
                    elif pos != last_pos + 1:
                        # Inside a valid region, but we've skipped a position.  Assume this means
                        # that the intervening positions are not covered at all.  Consequently,
                        # stay within state 1, but emit a line to denote the gap.
                        outfile.write('{}\n'.format(last_pos))
                        outfile.write('{}\t{}\t'.format(chrom, pos - 1))                        
                elif state == 2:
                    # Outside of a valid region. Verify that this is still the case
                    if frac >= minfrac:
                        # Entering a valid region.  Emit the start and transition to state 1.
                        outfile.write('{}\t{}\t'.format(chrom, pos - 1))
                        state = 1

                last_pos = pos

            if state == 1:
                outfile.write('{}\n'.format(last_pos + 1))

