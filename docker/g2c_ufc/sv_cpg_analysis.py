
def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--metadata', required=True, help='x_metadata.tsv')
    parser.add_argument('--rare_sv_tsv', required=True, help='sv phenotypes')
    parser.add_argument('--sample-list', required=True, help='list of samples to keep compared to larger G2C study')

if __name__ == '__main__':
    main()