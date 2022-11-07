import requests
import pandas as pd
import argparse
import random
import requests
import re
from bs4 import BeautifulSoup
import time


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


def get_record(pmcid):
    payload = {
        "verb": "GetRecord",
        "identifier": f"oai:pubmedcentral.nih.gov:{pmcid}",
        "metadataPrefix": "pmc",
    }
    r = requests.get("https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi", params=payload)
    return r


def pmc_query(ids, verbose=False):
    df = pd.DataFrame()
    for id in ids:
        try:
            rec = get_record(id)
        except TypeError as e:
            print("Could not find full text link: ", e)
            continue
        res = parse_record(rec.text)
        if verbose:
            print(id, res["detool"])
        df = pd.concat([df, pd.DataFrame.from_dict({id: res}, orient="index")])
    return df


def parse_record(xml):
    ex = re.compile(
        "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth"
    )
    soup = BeautifulSoup(xml, features="xml")
    if soup.find("error"):
        return {"detool": "", "context": soup.find("error").getText()}
    text = [p.getText() for p in soup.find_all("p")]
    matches = [ex.search(s.lower()) for s in text]
    detool = ";".join(set([m.group(0) for m in matches if m is not None]))
    context = "...".join([t for m, t in zip(matches, text) if m is not None])
    return {"detool": detool, "context": context}


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Scrapes differential expression analysis platform from Pubmed full text."
    )
    parser.add_argument("--infile", "-i", help="CSV file with PMCIDs", required=True)
    parser.add_argument(
        "--var",
        help="Variable name for PMCIDs. Defaults to PMCID.",
        default="PMCID",
    )
    parser.add_argument("--outfile", "-o", help="Output file", required=True)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-s", "--sample", action="store_true")
    parser.add_argument("--sleep", action="store_true")

    args = parser.parse_args()

    df = pd.read_csv(args.infile)
    ids = df[args.var].tolist()
    ids = [i for i in ids if str(i) != "nan"]
    ids = list(set(ids))
    if args.sample:
        ids = random.sample(ids, 5)

    chunks = list(divide_chunks(ids, 100))

    for count, chunk in enumerate(chunks):
        df = pmc_query(chunk, verbose=args.verbose)
        df.to_csv(args.outfile + str(count) + ".csv", index_label="PMCID")
        if args.sleep:
            nap = random.sample(range(60, 300, 1), 1)[0]
            print(f"Sleeping {nap}s, hold on..")
            time.sleep(nap)
