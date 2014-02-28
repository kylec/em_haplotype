import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-o' , dest='output')
args = parser.parse_args()


# filter input
def main():
    fin = open(args.input, 'r')
    fout = open(args.output, 'w')
    fout.write(fin.readline())
    
    for line in fin:
        fields = line.rstrip('\n').split('\t')
        keep = 1
        for field in fields:
            # remove individual with missing genotypes
            if (field == '.') | (field == ''):
                keep = 0
        if keep:
            fout.write(line)
                
    fin.close()
    fout.close()

if __name__ == '__main__':
    main()