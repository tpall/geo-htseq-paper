from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from bs4 import BeautifulSoup
import re
import pandas as pd
import argparse
import random
import time


class Pubmed:
    def __init__(self):
        self.chrome_options = Options()
        self.chrome_options.add_argument("--headless")
        self.driver = webdriver.Chrome(
            ChromeDriverManager().install(), options=self.chrome_options
        )
        self.detool = re.compile(
            "deseq2?|de(g|x)seq|rockhopper|cuff(diff|links)|edger|clc(bio)? ?genomics|igeak|bayseq|samseq|noiseq|bayseq|ebseq|limma|voom|sleuth|partek|(nrsa|nascent rna seq)|median ratio norm|rmats|ballgown|biojupie|seurat|exdega"
        )

    def close(self):
        self.driver.close()

    def get_detool(self, id):
        self.driver.get("https://doi.org/%s" % id)
        soup = BeautifulSoup(self.driver.page_source, "html.parser")
        if soup.find("error"):
            return {id: {"detool": "", "context": soup.find("error").getText()}}
        text = [p.getText() for p in soup.find_all("p")]
        matches = [
            ";".join(set([i.group() for i in self.detool.finditer(s.lower())]))
            for s in text
        ]
        detool = ";".join(
            set(";".join([m for m in matches if m is not None]).split(";"))
        )
        context = "...".join([t for m, t in zip(matches, text) if m is not None])
        return {id: {"detool": re.sub("^;", "", detool), "context": context}}

    def detool_query(self, ids, verbose=False):
        df = pd.DataFrame()
        for count, id in enumerate(ids):
            if count > 0:
                time.sleep(random.sample(range(1, 10, 1), 1)[0])
            try:
                res = self.get_detool(id)
            except TypeError as e:
                print("Could not find full text link: ", e)
                continue
            if verbose:
                print(id, res[id]["detool"])
            df = pd.concat([df, pd.DataFrame.from_dict(res, orient="index")])
        return df


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Scrapes differential expression analysis platform from Pubmed full text."
    )
    parser.add_argument("--infile", "-i", help="CSV file with DOIs", required=True)
    parser.add_argument(
        "--var",
        help="Variable name for DOI. Defaults to DOI.",
        default="DOI",
    )
    parser.add_argument("--outfile", "-o", help="Output file", required=True)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-s", "--sample", action="store_true")

    args = parser.parse_args()

    df = pd.read_csv(args.infile)
    ids = df[args.var].tolist()
    ids = [i for i in ids if str(i) != "nan"]
    ids = list(set(ids))
    if args.sample:
        ids = random.sample(ids, 5)

    chunks = list(divide_chunks(ids, 100))

    pm = Pubmed()
    for count, chunk in enumerate(chunks):
        if count > 0:
            nap = random.sample(range(240, 480, 1), 1)[0]
            print(f"Sleeping {nap}s, hold on..")
            time.sleep(nap)
        df = pm.detool_query(chunk, verbose=args.verbose)
        df.to_csv(args.outfile + str(count) + ".csv", index_label="PubMedIds")

    pm.close()
