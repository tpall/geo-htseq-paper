import requests
import pandas as pd
import argparse
import random
import requests


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


def pmid_to_pmcid(pmid, email, format):
    payload = {
        "tool": "python",
        "email": email,
        "ids": pmid,
        "idtype": "pmid",
        "format": format,
    }
    r = requests.get(
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/", params=payload
    )
    return r


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Scrapes differential expression analysis platform from Pubmed full text."
    )
    parser.add_argument("--infile", "-i", help="CSV file with PubMedIds", required=True)
    parser.add_argument(
        "--var",
        help="Variable name for PubMedIds. Defaults to PubMedIds.",
        default="PubMedIds",
    )
    parser.add_argument("--outfile", "-o", help="Output file", required=True)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-s", "--sample", action="store_true")
    parser.add_argument("-e", "--email", help="Your email address", required=True)
    parser.add_argument(
        "-f", "--format", help="Output format, defaults to csv", default="csv"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.infile)
    ids = df[args.var].tolist()
    ids = [i for i in ids if str(i) != "nan"]
    ids = list(set(ids))
    if args.sample:
        ids = random.sample(ids, 5)

    chunks = list(divide_chunks(ids, 200))

    for count, chunk in enumerate(chunks):
        pmids = ",".join(str(e) for e in chunk)
        if args.verbose:
            print("Working on chunk ", count)
        r = pmid_to_pmcid(pmid=pmids, email=args.email, format=args.format)
        with open(args.outfile + str(count) + ".csv", "w") as f:
            f.write(r.text)
