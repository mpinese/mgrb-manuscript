#!/usr/bin/env python3

if __name__ == '__main__':
    import sys

    sys.stdout.write(sys.stdin.readline())
    for line in sys.stdin:
        lineparts = line.rstrip().split('\t', 1)
        sys.stdout.write(lineparts[0])
        for dosage in lineparts[1]:
            sys.stdout.write('\t')
            if dosage == '.':
                sys.stdout.write('NA')
            else:
                sys.stdout.write(dosage)
        sys.stdout.write('\n')

