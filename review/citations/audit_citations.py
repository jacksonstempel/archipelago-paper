"""Citation audit for the ArchipelaGO language-model label review (test + development).

Stages (each cached; rerun is cheap):
  extract  -> citations.jsonl, decisions.jsonl   (from the two resolved_reviews.json files)
  resolve  -> resolves every identifier against public APIs (NCBI E-utilities, PMC idconv,
              doi.org handle API, Crossref, Europe PMC, UniProt REST, RCSB GraphQL, plain HTTP GET)
  match    -> matches each resolved record's text against the decision's protein / gene / organism
  report   -> summary.json, summary.md, failures.tsv, citation_results.tsv.gz

Usage: python -I audit_citations.py (--release DIR | --data DIR) --cache DIR --out DIR [--limit N] stage [stage ...]

--release DIR: the downloaded ArchipelaGO dataset (v2.0); the decisions are the review rows of annotations/{development,test}.parquet
and each decision's protein accessions those of proteins/{development,test}.parquet. --data DIR: the review records as the review wrote them
(DIR/test_resolved_reviews.json and DIR/dev_resolved_reviews.json, checked by their SHA-256), which the paper's audit read;
they are not part of the release. The released tables keep, for each decision, the URLs of the sources of the record that
settled it; the review records also hold the identifiers written in each source's reference text, the supplementary and
accessible URLs, the sources of the first review where an adjudication settled the decision, and the reviewer's stated
studied protein, title and year of each source (for the secondary checks). Run from the release, the audit therefore sees a
subset of the citations, and the secondary checks have nothing to compare. Matching uses the decision's accessions:
from the release, every UniProt accession of the protein's sequence; from the review records, those the records name.
Raw API responses are cached under --cache; an uncached run takes tens of minutes (public APIs, rate-limited).

Supplementary Table 8 is this audit of the review records (6 November 2026): 6,360 decisions, 6,358 citing, 11,548
citations of 6,872 distinct identifiers; 97.2% resolve, 80.6% of the resolving citations with text name the protein, 87.3% of
the citing decisions have such a citation. Run from the release (dataset v1.1, with the same API cache, 10 November 2026):
6,360 decisions, 6,358 citing, 10,735 citations of 6,358 distinct identifiers; 97.0% resolve, 81.7% name the protein, 87.6%
of the citing decisions have such a citation.
"""
import argparse, collections, concurrent.futures as cf, gzip, hashlib, html, json, os, random, re
import sys, threading, time, urllib.error, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

UA = "ArchipelaGO-citation-audit/1.0 (python-urllib)"
EXPECTED_SHA = {
    "test": "334bdb3115002bcd16c0e7ab0fd3010fada92f35c6c8354e4ed12210b8c017f1",
    "dev": "416cfbf34ac33830ca0b743edb3276da6720653f7a04703d71045e333c1003cf",
}
ARGS = None


# ----------------------------------------------------------------------------- utilities
def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


class Progress:
    def __init__(self, name, total):
        self.name, self.total, self.done, self.t0 = name, total, 0, time.time()
        self.last = 0

    def step(self, n=1):
        self.done += n
        now = time.time()
        if now - self.last > 10 or self.done == self.total:
            self.last = now
            el = now - self.t0
            rate = self.done / el if el > 0 else 0
            eta = (self.total - self.done) / rate if rate > 0 else float("nan")
            log(f"[{self.name}] {self.done}/{self.total} {rate:.1f}/s ETA {eta/60:.1f} min")


class RateLimiter:
    def __init__(self, per_sec):
        self.gap, self.lock, self.next = 1.0 / per_sec, threading.Lock(), 0.0

    def wait(self):
        with self.lock:
            now = time.time()
            if now < self.next:
                time.sleep(self.next - now)
            self.next = max(now, self.next) + self.gap


LIMITS = {"ncbi": RateLimiter(2.5), "crossref": RateLimiter(5), "epmc": RateLimiter(5),
          "uniprot": RateLimiter(5), "rcsb": RateLimiter(3), "doiorg": RateLimiter(5)}


def cache_path(kind, key):
    h = hashlib.sha1(key.encode()).hexdigest()
    d = os.path.join(ARGS.cache, kind, h[:2])
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, h + ".json")


def http(url, data=None, limiter=None, timeout=40, retries=4, max_bytes=None, headers=None):
    """Returns dict(status, final_url, ctype, body(bytes) or None, error)."""
    hdr = {"User-Agent": UA, "Accept": "*/*"}
    if headers:
        hdr.update(headers)
    last = None
    for attempt in range(retries):
        if limiter:
            limiter.wait()
        req = urllib.request.Request(url, data=data, headers=hdr)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read(max_bytes) if max_bytes else r.read()
                return {"status": r.status, "final_url": r.geturl(),
                        "ctype": r.headers.get("Content-Type", ""), "body": body, "error": None}
        except urllib.error.HTTPError as e:
            try:
                body = e.read(max_bytes or 2_000_000)
            except Exception:
                body = b""
            last = {"status": e.code, "final_url": e.geturl() if hasattr(e, "geturl") else url,
                    "ctype": e.headers.get("Content-Type", "") if e.headers else "", "body": body,
                    "error": f"HTTP {e.code}"}
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1 and limiter is not None:
                time.sleep(2 * (attempt + 1))
                continue
            return last
        except Exception as e:  # timeouts, DNS, TLS
            last = {"status": None, "final_url": url, "ctype": "", "body": None,
                    "error": f"{type(e).__name__}: {e}"[:300]}
            if attempt < retries - 1 and limiter is not None:
                time.sleep(2 * (attempt + 1))
                continue
            return last
    return last


def cached_json(kind, key, fn):
    p = cache_path(kind, key)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    v = fn()
    with open(p, "w", encoding="utf-8") as f:
        json.dump(v, f)
    return v


# ----------------------------------------------------------------------------- extraction
ACC_RE = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>]+")
PMID_TXT_RE = re.compile(r"PMID[:\s]*?(\d{5,9})", re.I)
PMC_RE = re.compile(r"PMC\d{4,9}")
JOURNAL_HOSTS = ("nature.com", "oup.com", "sciencedirect.com", "plos.org", "elifesciences.org",
                 "springer.com", "wiley.com", "jbc.org", "biologists.com", "frontiersin.org",
                 "mdpi.com", "cshlp.org", "asm.org", "rupress.org", "acs.org", "pnas.org",
                 "jstage.jst.go.jp", "jlr.org", "jci.org", "biorxiv.org", "physiology.org",
                 "tandfonline.com", "haematologica.org", "aacrjournals.org", "microbiologyresearch.org",
                 "spandidos-publications.com", "science.org", "cell.com", "embopress.org", "sciencemag.org",
                 "biomedcentral.com", "pubs.rsc.org", "annualreviews.org", "karger.com", "portlandpress.com",
                 "jneurosci.org", "royalsocietypublishing.org", "elsevier.com", "pnas", "ashpublications.org",
                 "ahajournals.org", "aspetjournals.org", "jimmunol.org", "genetics.org", "academic.oup",
                 "nar.oxfordjournals.org", "sagepub.com", "medrxiv.org", "thieme-connect", "degruyter.com",
                 "researchsquare.com", "peerj.com", "hindawi.com", "cambridge.org", "scielo", "jstor.org")
DB_HOSTS = ("ebi.ac.uk", "expasy.org", "informatics.jax.org", "brenda-enzymes", "flybase.org", "kegg.jp",
            "qmul.ac.uk", "reactome.org", "cazy.org", "yeastgenome.org", "zfin.org", "pombase.org",
            "rhea-db.org", "thebiogrid.org", "pdbj.org", "geneontology.org", "smart.embl.de",
            "inparanoidb", "pubchem.ncbi", "biocyc", "ncbi.nlm.nih.gov/gene", "ncbi.nlm.nih.gov/protein",
            "footprintdb", "ibdc.dbt.gov.in", "wormbase.org", "arabidopsis.org", "rgd.mcw.edu",
            "ensembl.org", "string-db.org", "interpro", "pfam", "alphafold", "dictybase", "candidagenome",
            "aureowiki", "pftargetbrowser", "seq2fun", "chemsensim", "xenbase", "tair", "plasmodb",
            "tritrypdb", "ecocyc", "metacyc", "rcsb.org", "uniprot.org", "mgi")


def host_of(u):
    m = re.match(r"https?://([^/?#]+)", u.strip(), re.I)
    return m.group(1).lower() if m else ""


def clean_doi(d):
    d = urllib.parse.unquote(d).strip()
    d = re.split(r"[?#]", d)[0]
    d = re.split(r"/mediaobjects/|_reference|/tables/|/figures/", d, flags=re.I)[0]
    d = re.sub(r"(\d{4}\.\d{2}\.\d{2}\.\d+)v\d+$", r"\1", d)  # preprint version suffix (bioRxiv-style)
    d = re.sub(r"(/(full|abstract|pdf|epdf|fulltext|full\.pdf|meta|figures?|suppl.*|references|"
               r"article-info|metrics))+/?$", "", d, flags=re.I)
    d = re.sub(r"\.(full|abstract|short|pdf|full\.pdf|full-text|long)(\.pdf)?$", "", d, flags=re.I)
    if d.startswith("10.1101/"):  # bioRxiv version suffix
        d = re.sub(r"v\d+$", "", d)
    return d.rstrip(".,;)]}").lower()


def classify_url(u):
    """-> (id_type, id_value, host_class)"""
    u = u.strip()
    h = host_of(u)
    path = urllib.parse.unquote(u)
    m = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d{1,9})(?:[/?#]|$)", path) or \
        re.search(r"ncbi\.nlm\.nih\.gov/pubmed/(\d{1,9})(?:[/?#]|$)", path) or \
        re.search(r"europepmc\.org/(?:article|abstract)/MED/(\d{1,9})", path, re.I)
    if m:
        return "pmid", m.group(1), "pubmed"
    m = re.search(r"(?:ncbi\.nlm\.nih\.gov/pmc/articles|pmc\.ncbi\.nlm\.nih\.gov/articles)/(PMC\d+)", path, re.I) or \
        re.search(r"europepmc\.org/(?:article|abstract)/PMC/(PMC\d+)", path, re.I) or \
        re.search(r"europepmc\.org/articles/(PMC\d+)", path, re.I)
    if m:
        return "pmcid", m.group(1).upper(), "pmc"
    if re.search(r"(dx\.)?doi\.org/", h + "/"):
        m = DOI_RE.search(path)
        if m:
            return "doi", clean_doi(m.group(0)), "doi"
    # Modern nature.com slugs (s41467-019-12345-6) are the DOI suffix; older slugs (nn0898_271, onc2013519,
    # 1202475) are not DOIs, so those stay plain URLs.
    m = re.search(r"nature\.com/articles/(s\d{5}-\d{3}-\d{4,5}-[0-9a-z])", path, re.I)
    if m:
        return "doi", clean_doi("10.1038/" + m.group(1)), "journal"
    m = re.search(r"elifesciences\.org/articles/(\d+)", path)
    if m:
        return "doi", "10.7554/elife." + m.group(1), "journal"
    m = re.search(r"(?:www\.)?uniprot\.org/(?:uniprotkb|uniprot)/(" + ACC_RE + r")(?:[./?#]|$)", path) or \
        re.search(r"rest\.uniprot\.org/uniprotkb/(" + ACC_RE + r")(?:[./?#]|$)", path)
    if m:
        return "uniprot", m.group(1), "uniprot"
    m = re.search(r"rcsb\.org/(?:structure|pdb/explore\.do\?structureId=|3d-view)/?([0-9][A-Za-z0-9]{3})(?:[/?#&]|$)", path, re.I) or \
        re.search(r"ebi\.ac\.uk/pdbe/entry/pdb/([0-9][A-Za-z0-9]{3})(?:[/?#]|$)", path, re.I) or \
        re.search(r"pdbj\.org/mine/summary/([0-9][A-Za-z0-9]{3})", path, re.I)
    if m:
        return "pdb", m.group(1).upper(), "pdb"
    m = DOI_RE.search(path)
    if m and any(j in h for j in JOURNAL_HOSTS + ("doi",)) or (m and "/doi/" in path):
        return "doi", clean_doi(m.group(0)), "journal"
    if any(j in h for j in JOURNAL_HOSTS):
        return "url", u, "journal"
    if any(dbh in (h + urllib.parse.urlparse(u).path) for dbh in DB_HOSTS):
        return "url", u, "other_database"
    return "url", u, "other"


def category(r):
    if r["decision"] == "accept" and r["review_support_type"] == "experimental":
        return "accepted_experimental"
    if r["decision"] == "accept" and r["review_support_type"] == "inference":
        return "accepted_inference"
    if r["decision"] == "do_not_accept":
        return "rejected"
    return f"other:{r['decision']}/{r['review_support_type']}"


def release_records(split):
    """The released decisions of a split as review records: the fields stage_extract reads (decision, support type,
    accessions, the deciding record's sources as URLs). Dataset v2.0: the review rows of annotations/<split> (the
    literature review of the GOA assertions; the development set's review of model false positives has no decision), the
    accessions from proteins/<split>, the term names from the shipped ontology."""
    import gzip
    import pyarrow.parquet as pq
    name = {"test": "test", "dev": "development"}[split]
    accessions = {r["protein_id"]: r["accessions"] for r in
                  pq.read_table(os.path.join(ARGS.release, "proteins", f"{name}.parquet"),
                                columns=["protein_id", "accessions"]).to_pylist()}
    names, term = {}, None
    with gzip.open(os.path.join(ARGS.release, "ontology", "go-basic.obo.gz"), "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("id: "):
                term = line[4:].strip()
            elif line.startswith("name: ") and term:
                names.setdefault(term, line[6:].strip())
            elif line.startswith("["):
                term = None
    out = []
    for r in pq.read_table(os.path.join(ARGS.release, "annotations", f"{name}.parquet")).to_pylist():
        if r["decision"] is None:
            continue
        out.append({"assertion_id": f"{r['protein_id']}:{r['go_id']}", "go_id": r["go_id"], "decision": r["decision"],
                    "review_support_type": r["support_type"],
                    "pro_review": {"go_name": names[r["go_id"]], "accessions": list(accessions[r["protein_id"]]),
                                   "sources": [{"url": u} for u in r["source_urls"] or []]}})
    return out

def stage_extract():
    decisions, cites, notes = [], [], collections.Counter()
    for split in ("test", "dev"):
        if ARGS.release:
            recs = release_records(split)
        else:
            fn = os.path.join(ARGS.data, f"{split}_resolved_reviews.json")
            raw = open(fn, "rb").read()
            sha = hashlib.sha256(raw).hexdigest()
            if sha != EXPECTED_SHA[split]:
                sys.exit(f"{fn}: sha256 {sha} != the review records' {EXPECTED_SHA[split]}")
            recs = json.loads(raw)
        for i, r in enumerate(recs):
            p = r.get("pro_review") or {}
            accs = []
            for k in ("accessions", "protein_accessions"):
                for a in p.get(k) or []:
                    if a not in accs:
                        accs.append(a)
            for g in r.get("original_goa_support") or []:
                a = g.get("accession")
                if a and a not in accs:
                    accs.append(a)
            if not accs:
                notes[f"{split}: decision without any accession"] += 1
            did = f"{split}:{r['assertion_id']}"
            decisions.append({"decision_id": did, "split": split, "assertion_id": r["assertion_id"],
                              "category": category(r), "go_id": r["go_id"], "go_name": p.get("go_name"),
                              "accessions": accs, "times_reviewed": r.get("times_reviewed"),
                              "adjudicated": bool(r.get("adjudication"))})
            srcs = [("pro_review", s) for s in (p.get("sources") or [])]
            adj = r.get("adjudication")
            if adj:
                srcs += [("adjudication", s) for s in (adj.get("sources") or [])]
            seen = set()
            for si, (origin, s) in enumerate(srcs):
                if isinstance(s, str):
                    s = {"url": s}
                elif not isinstance(s, dict):
                    notes[f"{split}: source of type {type(s).__name__}"] += 1
                    continue
                studied = " | ".join(str(s[k]) for k in ("studied_protein_or_organism", "studied_protein",
                                     "studied_entity", "studied_system", "studied_protein_or_system")
                                     if s.get(k))
                claimed_title = s.get("title") if isinstance(s.get("title"), str) else None
                cyear = None
                for k in ("publication_year", "publication_year_or_date", "publication_date"):
                    if s.get(k):
                        m = re.search(r"(19|20)\d\d", str(s[k]))
                        if m:
                            cyear = int(m.group(0))
                            break
                found = []
                for k in ("url", "supplementary_url", "accessible_url"):
                    v = s.get(k)
                    if isinstance(v, str) and v.strip():
                        if not re.match(r"https?://", v.strip(), re.I):
                            # non-URL text in a url field: try identifiers inside it
                            found += [(f"{k}(text)",) + t for t in ids_in_text(v)]
                            if not ids_in_text(v):
                                notes[f"{split}: non-URL text in source {k}"] += 1
                                found.append((f"{k}(text)", "unparsed", v.strip(), "unparsed"))
                        else:
                            found.append((k,) + classify_url(v))
                    elif v is not None and not isinstance(v, str):
                        notes[f"{split}: non-string {k}"] += 1
                if not s.get("url"):
                    notes[f"{split}: source without url field"] += 1
                for k in ("reference", "citation", "source_id", "publication"):
                    v = s.get(k)
                    if isinstance(v, str):
                        for t in ids_in_text(v):
                            found.append((k,) + t)
                for field, typ, val, hc in found:
                    key = (typ, val)
                    if key in seen:
                        continue
                    seen.add(key)
                    cites.append({"decision_id": did, "split": split, "category": decisions[-1]["category"],
                                  "origin": origin, "source_index": si, "field": field, "id_type": typ,
                                  "id": val, "host_class": hc, "host": host_of(val) if typ == "url" else "",
                                  "raw_url": s.get("url") if isinstance(s.get("url"), str) else None,
                                  "claimed_title": claimed_title, "claimed_year": cyear, "studied": studied,
                                  "section": s.get("section") or s.get("locator") or s.get("location") or
                                  s.get("section_or_record") or s.get("relevant_section_or_record") or ""})
    with open(os.path.join(ARGS.out_tmp, "decisions.jsonl"), "w", encoding="utf-8") as f:
        for d in decisions:
            f.write(json.dumps(d) + "\n")
    with open(os.path.join(ARGS.out_tmp, "citations.jsonl"), "w", encoding="utf-8") as f:
        for c in cites:
            f.write(json.dumps(c) + "\n")
    uniq = collections.Counter((c["id_type"]) for c in {(c["id_type"], c["id"]): c for c in cites}.values())
    log("decisions", len(decisions), "citation instances", len(cites), "unique ids", dict(uniq))
    log("notes", dict(notes))
    json.dump(dict(notes), open(os.path.join(ARGS.out_tmp, "extract_notes.json"), "w"), indent=1)


def ids_in_text(v):
    out = []
    for m in PMID_TXT_RE.finditer(v):
        out.append(("pmid", m.group(1), "pubmed"))
    for m in PMC_RE.finditer(v):
        out.append(("pmcid", m.group(0).upper(), "pmc"))
    for m in DOI_RE.finditer(v):
        out.append(("doi", clean_doi(m.group(0)), "doi"))
    return out


def load_jsonl(name):
    with open(os.path.join(ARGS.out_tmp, name), encoding="utf-8") as f:
        return [json.loads(l) for l in f]


# ----------------------------------------------------------------------------- resolvers
def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def txt(el):
    return "".join(el.itertext()).strip() if el is not None else ""


def parse_pubmed_xml(body):
    out = {}
    root = ET.fromstring(body)
    for art in root.findall("PubmedArticle"):
        pmid = txt(art.find("MedlineCitation/PMID"))
        a = art.find("MedlineCitation/Article")
        title = txt(a.find("ArticleTitle")) if a is not None else ""
        abstract = " ".join(txt(x) for x in a.findall("Abstract/AbstractText")) if a is not None else ""
        mesh = [txt(x) for x in art.findall("MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName")]
        chem = [txt(x) for x in art.findall("MedlineCitation/ChemicalList/Chemical/NameOfSubstance")]
        kw = [txt(x) for x in art.findall("MedlineCitation/KeywordList/Keyword")]
        genes = [txt(x) for x in art.findall("MedlineCitation/GeneSymbolList/GeneSymbol")]
        year = txt(art.find("MedlineCitation/Article/Journal/JournalIssue/PubDate/Year")) or \
            txt(art.find("MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate"))[:4]
        ids = {x.get("IdType"): txt(x) for x in art.findall("PubmedData/ArticleIdList/ArticleId")}
        out[pmid] = {"title": title, "abstract": abstract, "extra": " ; ".join(mesh + chem + kw + genes),
                     "year": int(year) if year[:4].isdigit() else None, "doi": ids.get("doi"),
                     "pmcid": ids.get("pmc")}
    for art in root.findall("PubmedBookArticle"):
        pmid = txt(art.find("BookDocument/PMID"))
        title = txt(art.find("BookDocument/ArticleTitle")) or txt(art.find("BookDocument/Book/BookTitle"))
        abstract = " ".join(txt(x) for x in art.findall("BookDocument/Abstract/AbstractText"))
        out[pmid] = {"title": title, "abstract": abstract, "extra": "", "year": None, "doi": None, "pmcid": None}
    return out


def resolve_pmids(pmids):
    res = {}
    todo = sorted(set(pmids))
    prog = Progress("pubmed", len(todo))
    for ch in chunks(todo, 200):
        key = "efetch:" + ",".join(ch)

        def fn():
            r = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                     data=urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ch), "retmode": "xml",
                                                  "tool": "ArchipelaGO-citation-audit"}).encode(),
                     limiter=LIMITS["ncbi"], timeout=120)
            if r["status"] != 200:
                raise RuntimeError(f"efetch failed {r['status']} {r['error']}")
            return {"xml": r["body"].decode("utf-8", "replace")}
        got = parse_pubmed_xml(cached_json("efetch", key, fn)["xml"].encode("utf-8"))
        for p in ch:
            if p in got:
                res[p] = {"resolves": "yes", "how": "PubMed efetch", **got[p]}
            else:
                res[p] = {"resolves": "no", "how": "PubMed efetch returned no record"}
        prog.step(len(ch))
    # confirm the misses one more time with esummary (catches efetch hiccups)
    miss = [p for p in todo if res[p]["resolves"] == "no"]
    for ch in chunks(miss, 200):
        def fn():
            r = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "pubmed", "id": ",".join(ch), "retmode": "json"}), limiter=LIMITS["ncbi"])
            return json.loads(r["body"]) if r["status"] == 200 else {"error": r["error"]}
        js = cached_json("esummary", "pubmed:" + ",".join(ch), fn)
        for p in ch:
            e = (js.get("result") or {}).get(p)
            if e and not e.get("error"):
                res[p] = {"resolves": "yes", "how": "PubMed esummary", "title": e.get("title", ""),
                          "abstract": "", "extra": "", "year": None, "doi": None, "pmcid": None}
            elif e:
                res[p]["how"] += f"; esummary: {e.get('error')}"
    return res


def epmc_search(query):
    def fn():
        r = http("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
            {"query": query, "resultType": "core", "format": "json", "pageSize": 5}), limiter=LIMITS["epmc"])
        if r["status"] != 200:
            return {"error": r["error"], "status": r["status"]}
        return json.loads(r["body"])
    return cached_json("epmc", query, fn)


def epmc_record(js):
    lst = ((js.get("resultList") or {}).get("result")) or []
    if not lst:
        return None
    e = lst[0]
    mesh = [m.get("descriptorName", "") for m in ((e.get("meshHeadingList") or {}).get("meshHeading") or [])]
    kw = (e.get("keywordList") or {}).get("keyword") or []
    chem = [c.get("name", "") for c in ((e.get("chemicalList") or {}).get("chemical") or [])]
    y = e.get("pubYear")
    return {"title": e.get("title", "") or "", "abstract": re.sub(r"<[^>]+>", " ", e.get("abstractText", "") or ""),
            "extra": " ; ".join(mesh + kw + chem), "year": int(y) if y and y.isdigit() else None,
            "pmid": str(e["pmid"]) if e.get("pmid") else None, "pmcid": e.get("pmcid"), "doi": (e.get("doi") or "").lower() or None}


def resolve_pmcids(pmcids):
    res, todo = {}, sorted(set(pmcids))
    prog = Progress("pmc-idconv", len(todo))
    for ch in chunks(todo, 200):
        def fn():
            r = http("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?" + urllib.parse.urlencode(
                {"ids": ",".join(ch), "format": "json", "tool": "ArchipelaGO-citation-audit"}),
                limiter=LIMITS["ncbi"], timeout=90)
            if r["status"] != 200:
                raise RuntimeError(f"idconv {r['status']} {r['error']}")
            return json.loads(r["body"])
        js = cached_json("idconv", ",".join(ch), fn)
        for rec in js.get("records", []):
            q = (rec.get("requested-id") or rec.get("pmcid") or "").upper()
            if rec.get("status") == "error" or not rec.get("pmcid"):
                res[q] = {"resolves": "pending", "how": f"idconv: {rec.get('errmsg', 'error')}"}
            else:
                res[q] = {"resolves": "yes", "how": "PMC idconv",
                          "pmid": str(rec["pmid"]) if rec.get("pmid") else None,
                          "doi": (rec.get("doi") or "").lower() or None}
        prog.step(len(ch))
    # idconv misses: confirm with esummary db=pmc, then Europe PMC
    for p in todo:
        if p not in res or res[p]["resolves"] == "pending":
            def fn():
                r = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
                    {"db": "pmc", "id": p[3:], "retmode": "json"}), limiter=LIMITS["ncbi"])
                return json.loads(r["body"]) if r["status"] == 200 else {"error": r["error"]}
            js = cached_json("esummary", "pmc:" + p, fn)
            e = (js.get("result") or {}).get(p[3:])
            prev = (res.get(p) or {}).get("how", "idconv: not returned")
            if e and not e.get("error") and e.get("title"):
                ids = {a.get("idtype"): a.get("value") for a in e.get("articleids", [])}
                res[p] = {"resolves": "yes", "how": prev + "; PMC esummary", "title": e.get("title", ""),
                          "pmid": ids.get("pmid") if ids.get("pmid") not in (None, "0") else None,
                          "doi": (ids.get("doi") or "").lower() or None}
            else:
                er = epmc_record(epmc_search(f"PMCID:{p}"))
                if er:
                    res[p] = {"resolves": "yes", "how": prev + "; Europe PMC", **er}
                else:
                    res[p] = {"resolves": "no", "how": prev + f"; PMC esummary: {(e or {}).get('error', 'none')}; Europe PMC: none"}
    return res


def resolve_dois(dois):
    res, todo = {}, sorted(set(dois))
    prog = Progress("doi", len(todo))

    def one(d):
        def hfn():
            r = http("https://doi.org/api/handles/" + urllib.parse.quote(d, safe="/"), limiter=LIMITS["doiorg"])
            try:
                return {"status": r["status"], "js": json.loads(r["body"] or b"{}")}
            except Exception:
                return {"status": r["status"], "js": {}, "error": r["error"]}
        h = cached_json("doihandle", d, hfn)
        code = (h.get("js") or {}).get("responseCode")
        if code == 1:
            out = {"resolves": "yes", "how": "doi.org handle"}
        elif code == 100 or h.get("status") == 404:
            out = {"resolves": "no", "how": "doi.org handle: DOI not found"}
        else:
            out = {"resolves": "unknown", "how": f"doi.org handle status {h.get('status')} code {code}"}

        def cfn():
            r = http("https://api.crossref.org/works/" + urllib.parse.quote(d, safe="/"), limiter=LIMITS["crossref"])
            if r["status"] == 200:
                return json.loads(r["body"])
            return {"status": r["status"]}
        cr = cached_json("crossref", d, cfn)
        msg = cr.get("message") if isinstance(cr, dict) else None
        if msg:
            out["title"] = " ".join(msg.get("title") or [])
            out["abstract"] = re.sub(r"<[^>]+>", " ", msg.get("abstract") or "")
            out["extra"] = " ; ".join((msg.get("subject") or []) + (msg.get("subtitle") or []))
            dp = (msg.get("issued") or {}).get("date-parts") or [[None]]
            out["year"] = dp[0][0] if dp and dp[0] else None
            out["how"] += "; Crossref"
        er = epmc_record(epmc_search(f'DOI:"{d}"'))
        if er:
            for k in ("title", "abstract", "extra", "year"):
                if er.get(k) and (len(str(er[k])) > len(str(out.get(k) or ""))):
                    out[k] = er[k]
            out["pmid"] = er.get("pmid")
            out["how"] += "; Europe PMC"
        if out["resolves"] != "yes" and (msg or er):
            out["resolves"] = "yes"
        return d, out

    with cf.ThreadPoolExecutor(4) as ex:
        for d, out in ex.map(one, todo):
            res[d] = out
            prog.step()
    return res


def uniprot_fetch(accs):
    """accession -> dict(active, genes[], protein_names[], organism{sci, common, syn}, taxid, raw_ok)"""
    res, todo = {}, sorted(set(accs))
    prog = Progress("uniprot", len(todo))
    fields = "accession,id,protein_name,gene_names,organism_name,organism_id"
    for ch in chunks(todo, 200):
        def fn():
            r = http("https://rest.uniprot.org/uniprotkb/accessions?" + urllib.parse.urlencode(
                {"accessions": ",".join(ch), "format": "json", "fields": fields, "size": 500}),
                limiter=LIMITS["uniprot"], timeout=120)
            if r["status"] != 200:
                return {"error": r["error"], "status": r["status"],
                        "body": (r["body"] or b"")[:2000].decode("utf-8", "replace")}
            return json.loads(r["body"])
        js = cached_json("uniprot", ",".join(ch), fn)
        if "results" not in js:
            log("UniProt batch error", js.get("status"), js.get("body", "")[:300])
            # fall back to one-by-one so a single bad accession doesn't sink the batch
            for a in ch:
                js1 = cached_json("uniprot1", a, lambda a=a: (lambda r: json.loads(r["body"]) if r["status"] == 200 else
                                  {"error": r["error"], "status": r["status"]})(
                    http(f"https://rest.uniprot.org/uniprotkb/{a}.json?fields={fields}", limiter=LIMITS["uniprot"])))
                if "primaryAccession" in js1:
                    res[a] = parse_uniprot(js1)
                else:
                    res[a] = {"active": False, "status": js1.get("status"), "reason": js1.get("error")}
        else:
            for e in js["results"]:
                pe = parse_uniprot(e)
                res[e["primaryAccession"]] = pe
                for sa in e.get("secondaryAccessions") or []:
                    if sa in ch and sa not in res:
                        res[sa] = pe
            for a in ch:
                if a not in res:
                    # the batch endpoint omits some inactive (deleted/merged) entries; ask for the single entry
                    js1 = cached_json("uniprot1", a, lambda a=a: (lambda r: json.loads(r["body"]) if r["status"] == 200 else
                                      {"error": r["error"], "status": r["status"]})(
                        http(f"https://rest.uniprot.org/uniprotkb/{a}.json", limiter=LIMITS["uniprot"])))
                    if "primaryAccession" in js1:
                        res[a] = parse_uniprot(js1)
                        if not res[a].get("active"):
                            res[a]["reason"] = "inactive: " + res[a].get("reason", "")
                    else:
                        res[a] = {"active": False, "reason": f"not found in UniProt ({js1.get('error')})"}
        prog.step(len(ch))
    return res


def parse_uniprot(e):
    if e.get("entryType", "").startswith("Inactive"):
        return {"active": False, "reason": (e.get("inactiveReason") or {}).get("inactiveReasonType", "inactive"),
                "accession": e.get("primaryAccession")}
    genes = []
    for g in e.get("genes") or []:
        for k in ("geneName",):
            if g.get(k):
                genes.append(g[k]["value"])
        for k in ("synonyms", "orderedLocusNames", "orfNames"):
            genes += [x["value"] for x in g.get(k) or []]
    names = []
    pd = e.get("proteinDescription") or {}

    def add(n):
        if not n:
            return
        if n.get("fullName"):
            names.append(n["fullName"]["value"])
        names.extend(x["value"] for x in n.get("shortNames") or [])
    add(pd.get("recommendedName"))
    for n in pd.get("alternativeNames") or []:
        add(n)
    for n in pd.get("submissionNames") or []:
        add(n)
    for grp in ("includes", "contains"):
        for sub in pd.get(grp) or []:
            add(sub.get("recommendedName"))
            for n in sub.get("alternativeNames") or []:
                add(n)
    org = e.get("organism") or {}
    return {"active": True, "accession": e.get("primaryAccession"), "uniprot_id": e.get("uniProtkbId"),
            "genes": genes, "protein_names": names, "organism": org.get("scientificName", ""),
            "organism_common": org.get("commonName", ""), "organism_syn": org.get("synonyms") or [],
            "taxid": org.get("taxonId")}


def resolve_pdb(ids):
    res, todo = {}, sorted(set(ids))
    prog = Progress("pdb", len(todo))
    q = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title}
      polymer_entities{rcsb_polymer_entity{pdbx_description}
      rcsb_polymer_entity_container_identifiers{uniprot_ids}
      rcsb_entity_source_organism{scientific_name ncbi_common_names}
      rcsb_entity_host_organism{scientific_name}}}}"""
    for ch in chunks(todo, 50):
        def fn():
            r = http("https://data.rcsb.org/graphql", data=json.dumps({"query": q, "variables": {"ids": ch}}).encode(),
                     limiter=LIMITS["rcsb"], headers={"Content-Type": "application/json"}, timeout=90)
            if r["status"] != 200:
                raise RuntimeError(f"rcsb {r['status']} {r['error']}")
            return json.loads(r["body"])
        js = cached_json("rcsb", ",".join(ch), fn)
        got = {}
        for e in ((js.get("data") or {}).get("entries") or []):
            if not e:
                continue
            desc, unis, orgs = [], [], []
            for pe in e.get("polymer_entities") or []:
                desc.append(((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")) or "")
                unis += ((pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids")) or []
                for o in pe.get("rcsb_entity_source_organism") or []:
                    orgs.append(o.get("scientific_name") or "")
                    orgs += o.get("ncbi_common_names") or []
            got[e["rcsb_id"].upper()] = {"resolves": "yes", "how": "RCSB GraphQL",
                                         "title": (e.get("struct") or {}).get("title") or "",
                                         "abstract": " ; ".join(desc), "extra": " ; ".join(orgs),
                                         "uniprot_ids": sorted(set(unis))}
        for i in ch:
            res[i] = got.get(i) or {"resolves": "no", "how": "RCSB GraphQL: entry not found"}
        prog.step(len(ch))
    return res


def strip_html(b, ctype):
    s = b.decode("utf-8", "replace")
    m = re.search(r"<title[^>]*>(.*?)</title>", s, re.I | re.S)
    title = html.unescape(re.sub(r"\s+", " ", m.group(1))).strip() if m else ""
    metas = re.findall(r'<meta[^>]+(?:name|property)="(?:citation_title|dc\.title|og:title|description|'
                       r'og:description|citation_abstract|dc\.description)"[^>]+content="([^"]*)"', s, re.I)
    if metas and (not title or len(title) < 15):
        title = html.unescape(metas[0])
    s = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", s, flags=re.I | re.S)
    body = html.unescape(re.sub(r"<[^>]+>", " ", s))
    body = re.sub(r"\s+", " ", body)
    return title, " ".join(html.unescape(x) for x in metas) + " " + body


def pdf_text(b):
    try:
        import fitz
        doc = fitz.open(stream=b, filetype="pdf")
        t = " ".join(doc[i].get_text() for i in range(min(len(doc), 30)))
        title = (doc.metadata or {}).get("title") or ""
        return title, t
    except Exception as e:
        return "", ""


def ncbi_esummary_urls(urls):
    """NCBI Gene/Protein web pages answer scripted GETs with a reCAPTCHA page; use E-utilities instead."""
    res, want = {}, collections.defaultdict(dict)
    for u in urls:
        m = re.search(r"ncbi\.nlm\.nih\.gov/(gene|protein)/(?:\?term=)?([A-Za-z0-9_.]+)/?(?:[?#]|$)", u)
        if m:
            want[m.group(1)][u] = m.group(2)
    for db, um in want.items():
        ids = sorted(set(um.values()))
        got = {}
        for ch in chunks(ids, 200):
            def fn():
                r = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                         data=urllib.parse.urlencode({"db": db, "id": ",".join(ch), "retmode": "json"}).encode(),
                         limiter=LIMITS["ncbi"], timeout=90)
                return json.loads(r["body"]) if r["status"] == 200 else {"error": r["error"]}
            js = cached_json("esummary", f"{db}:" + ",".join(ch), fn)
            result = js.get("result") or {}
            for uid in result.get("uids", []):
                e = result[uid]
                keys = {uid, str(e.get("accessionversion", "")), str(e.get("caption", ""))}
                if db == "gene":
                    t = f"{e.get('name', '')} {e.get('description', '')}"
                    x = f"{e.get('otheraliases', '')} ; {e.get('otherdesignations', '')} ; " \
                        f"{(e.get('organism') or {}).get('scientificname', '')} ; {(e.get('organism') or {}).get('commonname', '')}"
                    ab = e.get("summary", "")
                else:
                    t, x, ab = e.get("title", ""), e.get("organism", ""), ""
                ok = not e.get("error")
                for k in keys:
                    got[k] = {"resolves": "yes" if ok else "no", "how": f"NCBI {db} esummary"
                              + (f" (status {e.get('status')}, current id {e.get('currentid')})" if e.get('currentid') else ""),
                              "title": t, "abstract": ab, "extra": x}
        for u, i in um.items():
            res[u] = got.get(i) or got.get(i.split(".")[0]) or {"resolves": "no", "how": f"NCBI {db} esummary: no record"}
    return res


def rewrite_for_fetch(u):
    """Single-page apps whose HTML holds no content: fetch the documented JSON endpoint instead."""
    m = re.search(r"ebi\.ac\.uk/interpro/(?:api/)?entry/(interpro|pfam|InterPro|Pfam|cathgene3d|panther|smart|"
                  r"prosite|ncbifam|cdd|hamap|pirsf|profile|sfld|ssf|prints)/([A-Za-z0-9_.]+)", u, re.I)
    if m:
        return f"https://www.ebi.ac.uk/interpro/api/entry/{m.group(1).lower()}/{m.group(2)}/"
    m = re.search(r"//www\.uniprot\.org/(taxonomy|uniref|uniparc|keywords|unirule|citations|diseases|locations)/([A-Za-z0-9_.-]+)", u)
    if m:
        return f"https://rest.uniprot.org/{m.group(1)}/{m.group(2)}.json"
    m = re.search(r"ebi\.ac\.uk/QuickGO/term/(GO:\d{7})", u)
    if m:
        return f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{m.group(1)}"
    return u


def resolve_urls(urls):
    res, todo = {}, sorted(set(urls))
    res.update(ncbi_esummary_urls(todo))
    todo = [u for u in todo if u not in res]
    by_host = collections.defaultdict(list)
    for u in todo:
        by_host[host_of(u)].append(u)
    prog = Progress("urls", len(todo))
    lock = threading.Lock()

    def fetch_host(h):
        out = {}
        for u in by_host[h]:
            fu = rewrite_for_fetch(u)

            def fn():
                r = http(fu, timeout=30, retries=2, max_bytes=8_000_000, limiter=RateLimiter(1.5),
                         headers={"Accept": "text/html,application/xhtml+xml,application/xml,application/json,"
                                            "application/pdf;q=0.9,*/*;q=0.8"})
                body, ctype = r["body"] or b"", r["ctype"] or ""
                title, text = "", ""
                if body:
                    if "pdf" in ctype or body[:5] == b"%PDF-":
                        title, text = pdf_text(body)
                    elif "json" in ctype:
                        text = body.decode("utf-8", "replace")
                    else:
                        title, text = strip_html(body, ctype)
                return {"status": r["status"], "final_url": r["final_url"], "ctype": ctype, "error": r["error"],
                        "title": title[:500], "text": text[:400_000]}
            r = cached_json("url", fu, fn)
            st = r["status"]
            if st and 200 <= st < 300:
                rv = "yes"
            elif st in (401, 403, 429, 451) or (st in (202,)):
                rv = "blocked"
            elif st in (404, 410):
                rv = "no"
            else:
                rv = "unreachable"
            # Soft 404 detection
            if rv == "yes" and re.search(r"\b(page not found|404 not found|article not found|no longer available|"
                                         r"we can.t find|cannot be found|does not exist)\b", (r["title"] or "") + " " +
                                         ((r["text"] or "") if len(r["text"] or "") < 3000 else ""), re.I):
                rv = "soft404?"
            # a deep link that redirects to a site's root or home page is a dead link
            if rv == "yes" and r["final_url"] and r["final_url"] != fu:
                fp = urllib.parse.urlparse(r["final_url"]).path.rstrip("/")
                if fp in ("", "/index.html", "/index.php", "/home") and len(urllib.parse.urlparse(fu).path) > 1:
                    rv = "soft404?"
            if rv == "yes" and re.search(r"(just a moment|attention required|captcha|access denied|are you a robot|"
                                         r"enable javascript and cookies)", (r["title"] or "") + " " + (r["text"] or "")[:2000], re.I):
                rv = "blocked"
            out[u] = {"resolves": rv, "how": (f"HTTP GET {st} {r['error'] or ''}".strip() +
                                              (f" (fetched {fu})" if fu != u else "")), "title": r["title"],
                      "abstract": "", "page_text": r["text"], "final_url": r["final_url"], "ctype": r["ctype"]}
            with lock:
                prog.step()
        return out
    with cf.ThreadPoolExecutor(12) as ex:
        for out in ex.map(fetch_host, sorted(by_host, key=lambda h: -len(by_host[h]))):
            res.update(out)
    # Elsevier-style URLs (sciencedirect /pii/, jbc.org, jlr.org, cell.com ...) embed a PII: look it up in Crossref
    # (alternative-id) so paywalled pages still get a title/abstract record.
    piis = {}
    for u in todo:
        m = re.search(r"(?:/pii/|/article/)(S?[0-9X()\-]{15,25}[0-9X])", urllib.parse.unquote(u), re.I)
        if m:
            piis[u] = re.sub(r"[^0-9A-Za-z]", "", m.group(1)).upper()
    if piis:
        prog = Progress("pii-crossref", len(piis))
        for u, pii in piis.items():
            def fn():
                r = http("https://api.crossref.org/works?" + urllib.parse.urlencode(
                    {"filter": f"alternative-id:{pii}", "rows": 2}), limiter=LIMITS["crossref"])
                return json.loads(r["body"]) if r["status"] == 200 else {"status": r["status"]}
            js = cached_json("crossref_pii", pii, fn)
            items = ((js.get("message") or {}).get("items")) or []
            v = res[u]
            v["pii"] = pii
            if items:
                d = (items[0].get("DOI") or "").lower()
                v["pii_doi"] = d
                dr = resolve_dois([d])[d]
                v["resolves"] = "yes"
                v["how"] += f"; PII {pii} -> Crossref DOI {d} ({dr['how']})"
                for k in ("title", "abstract", "extra", "year"):
                    if dr.get(k):
                        v[k] = dr[k]
            else:
                v["how"] += f"; PII {pii} not in Crossref"
            prog.step()
    return res


def stage_resolve():
    cites = load_jsonl("citations.jsonl")
    decs = load_jsonl("decisions.jsonl")
    if ARGS.limit:
        keep = set(d["decision_id"] for d in random.Random(0).sample(decs, ARGS.limit))
        cites = [c for c in cites if c["decision_id"] in keep]
        decs = [d for d in decs if d["decision_id"] in keep]
    by = collections.defaultdict(set)
    for c in cites:
        by[c["id_type"]].add(c["id"])
    log("unique ids:", {k: len(v) for k, v in by.items()})
    est = (len(by["pmid"]) / 200 + len(by["pmcid"]) / 200 * 2) / 2.5 + len(by["doi"]) * 3 / 5 / 4 * 4 + \
        len(by["uniprot"]) / 200 + len(by["pdb"]) / 50 + len(by["url"]) / 6
    log(f"rough runtime estimate (uncached): {est/60:.1f} min")
    R = {}
    R["pmcid"] = resolve_pmcids(by["pmcid"])
    # add PMIDs behind PMCIDs and DOIs for title/abstract
    extra_pmids = {v.get("pmid") for v in R["pmcid"].values() if v.get("pmid")}
    R["pdb"] = resolve_pdb(by["pdb"])
    R["doi"] = resolve_dois(by["doi"])
    extra_pmids |= {v.get("pmid") for v in R["doi"].values() if v.get("pmid")}
    R["pmid"] = resolve_pmids(by["pmid"] | {p for p in extra_pmids if p})
    for p, v in R["pmcid"].items():
        if v.get("pmid") and v["pmid"] in R["pmid"] and R["pmid"][v["pmid"]]["resolves"] == "yes":
            pm = R["pmid"][v["pmid"]]
            for k in ("title", "abstract", "extra", "year"):
                if pm.get(k) and not v.get(k):
                    v[k] = pm[k]
            v["how"] += "; PubMed efetch via PMID"
        elif v["resolves"] == "yes" and not v.get("abstract"):
            er = epmc_record(epmc_search(f"PMCID:{p}"))
            if er:
                for k in ("title", "abstract", "extra", "year"):
                    if er.get(k) and not v.get(k):
                        v[k] = er[k]
                v["how"] += "; Europe PMC"
    for d, v in R["doi"].items():
        if v.get("pmid") and v["pmid"] in R["pmid"] and R["pmid"][v["pmid"]]["resolves"] == "yes":
            pm = R["pmid"][v["pmid"]]
            for k in ("title", "abstract", "extra"):
                if pm.get(k) and len(pm[k]) > len(v.get(k) or ""):
                    v[k] = pm[k]
    # UniProt: decision accessions + cited accessions
    dec_accs = {a for d in decs for a in d["accessions"]}
    U = uniprot_fetch(dec_accs | by["uniprot"])
    R["uniprot"] = {}
    for a in by["uniprot"]:
        u = U.get(a) or {}
        if u.get("active"):
            R["uniprot"][a] = {"resolves": "yes", "how": "UniProt REST", "title": "; ".join(u["protein_names"][:3]),
                               "abstract": "", "extra": " ; ".join(u["genes"] + [u["organism"], u["organism_common"]]),
                               "uniprot": u}
        elif "accession" in u and not u.get("active"):   # entry exists but is inactive (deleted/merged)
            R["uniprot"][a] = {"resolves": "yes", "how": f"UniProt REST: inactive entry ({u.get('reason')})",
                               "title": "", "abstract": "", "extra": "", "inactive": True}
        else:
            R["uniprot"][a] = {"resolves": "no", "how": f"UniProt REST: {u.get('reason')}"}
    R["url"] = resolve_urls(by["url"])
    R["unparsed"] = {u: {"resolves": "no", "how": "not a URL or identifier"} for u in by["unparsed"]}
    with gzip.open(os.path.join(ARGS.out_tmp, "resolved.json.gz"), "wt", encoding="utf-8") as f:
        json.dump(R, f)
    with gzip.open(os.path.join(ARGS.out_tmp, "uniprot_decisions.json.gz"), "wt", encoding="utf-8") as f:
        json.dump(U, f)
    log("resolve done:", {k: collections.Counter(v["resolves"] for v in R[k].values()) for k in R})


# ----------------------------------------------------------------------------- matching
GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa", "λ": "lambda",
         "μ": "mu", "σ": "sigma", "ζ": "zeta", "η": "eta", "θ": "theta", "ω": "omega", "τ": "tau", "ρ": "rho"}
GENE_STOP = set("""the and for not set cat ace arc max ran was can has age all ask big bad but cut end far fat fit gap
hit ice key lab let map mix net out pan per pet pin put red sun tag top two via war yes dna rna atp adp amp gtp gdp nad
fad fmn coa pcr ros gfp his cell type wild mrna trna rrna camp cgmp nadh nadp kda one ten age aim air arm art bag bat bed
bee bin bit box boy bud bus cap car cop cow cry cup day den dim dog dot dry ear egg ego elf elk fan fig fly fog fox fun fur
gas gel gem gin god gum gun gut hat hay hen hip hog hop hot hub hug ion jam jar jaw jet job jog joy lag lap law lay leg lid
lie lip log low mad man mat men met mob mod mop mud mug nap net new nod nor now nut oak odd oil old orb ore owl own pad pal
paw pay pea pen pie pig pit pod pop pot pro pub pup ram rap rat raw ray rib rid rim rip rod rot row rub rug run sac sad saw
sea see sex shy sin sip sit six ski sky sly sob son soy spa spy sub sum tab tan tap tar tea tee tie tin tip toe ton toy try
tub tug use van vat vet vow wax way web wet who why wig win wit wok won yak yam yet zip zoo also than that this with from
into were been have which their there when what will more most some such only other these those then them they here very
loss gain long short light dark heavy rest wing eye eyes hair tail head body foot hand heart brain blood bone skin muscle
early late white black green yellow brown blue gold iron lead zinc copper silver salt sugar fat fats acid base virus
protein proteins gene genes factor kinase domain binding complex subunit receptor channel enzyme""".split())
GENERIC_WORDS = set("""protein proteins putative probable uncharacterized uncharacterised kinase factor subunit domain
containing family member like homolog homologue chain alpha beta gamma delta epsilon receptor transporter binding enzyme
isoform fragment precursor type class related associated regulator regulatory transcription transcriptional dna rna
small large component complex unit predicted hypothetical conserved unknown function activity dependent specific
non and of the a an in to with on for or by ii iii iv v vi i atp gtp nad nadp nadh coa cell surface membrane nuclear
mitochondrial cytoplasmic periplasmic outer inner integral transmembrane integral sensor response system two
serine threonine tyrosine protein-kinase phosphatase synthase synthetase reductase oxidase dehydrogenase transferase
hydrolase ligase lyase isomerase mutase esterase lipase protease peptidase ribosomal ribosome rrna trna mrna
heat shock zinc finger finger box helix loop leucine rich repeat repeats wd ankyrin motif region terminal n c
transport export import pump permease atpase gtpase methyltransferase acetyltransferase glycosyltransferase
long short chain cytochrome b c p450 polymerase helicase nuclease endonuclease exonuclease""".split())
COMMON_PLURAL = {"mouse": ["mice"], "human": ["humans"], "rat": ["rats"], "chicken": ["chickens"],
                 "zebrafish": [], "baker's yeast": ["yeast", "budding yeast"], "fission yeast": ["yeast"],
                 "fruit fly": ["drosophila", "flies"], "thale cress": ["arabidopsis"], "bovine": ["cattle", "cow"],
                 "cattle": ["bovine", "cows"], "pig": ["porcine", "swine", "pigs"], "dog": ["canine", "dogs"],
                 "rabbit": ["rabbits"], "sheep": ["ovine"], "horse": ["equine"],
                 "african clawed frog": ["xenopus"], "rice": ["oryza"], "maize": ["zea mays", "corn"]}


def norm(s):
    s = (s or "").lower()
    for k, v in GREEK.items():
        s = s.replace(k, v)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"(?<=[a-z]) (?=\d)", "", s)   # "pen 2" -> "pen2"
    s = re.sub(r"(?<=\d) (?=[a-z]\b)", "", s)  # "pen2 a" -> "pen2a" (single trailing letter)
    s = re.sub(r"(?<=\d) (alpha|beta|gamma|delta|epsilon|kappa)\b", r"\1", s)  # "pgc1 alpha" -> "pgc1alpha"
    return " " + re.sub(r"\s+", " ", s).strip() + " "


def phrase_in(p, ntext):
    """p already normalized (with surrounding spaces)."""
    return p in ntext


def gene_variants(g):
    n = norm(g).strip()
    if len(n.replace(" ", "")) < 3 or n in GENE_STOP or n.isdigit():
        return []
    return [" " + n + " ", " " + n + "p "]   # allow yeast 'p' protein suffix (Cdc42p)


def protein_variants(name):
    n = norm(name).strip()
    toks = n.split()
    if not toks:
        return []
    distinct = [t for t in toks if t not in GENERIC_WORDS and not t.isdigit() and len(t) >= 3]
    if not distinct:
        return []
    if len(toks) == 1 and (len(n) < 4 or n in GENE_STOP):
        return []
    out = [" " + n + " "]
    m = re.fullmatch(r"([a-z]{5,})\d{1,2}", n)          # "hepcidin1" -> also "hepcidin"
    if m and m.group(1) not in GENERIC_WORDS:
        out.append(" " + m.group(1) + " ")
    m = re.fullmatch(r"([a-z]{4,}) (alpha|beta|gamma|delta|epsilon|kappa)\d{0,2}", n)  # "catenin beta1" -> "beta catenin"
    if m and m.group(1) not in GENERIC_WORDS:
        out.append(f" {m.group(2)} {m.group(1)} ")
    return out


def organism_variants(u):
    out = set()
    sci = re.sub(r"\(.*?\)", "", u.get("organism") or "").strip()
    parts = sci.split()
    if len(parts) >= 2:
        out.add(" ".join(parts[:2]))
        out.add(parts[0][0] + " " + parts[1])         # "s cerevisiae" (from "S. cerevisiae")
        if len(parts[0]) >= 6:
            out.add(parts[0])                           # genus
    elif parts:
        out.add(parts[0])
    if sci:
        out.add(sci)
    com = (u.get("organism_common") or "").lower()
    for c in [com] + [x.lower() for x in u.get("organism_syn") or []]:
        if c and len(c) >= 3:
            out.add(c)
            out.update(COMMON_PLURAL.get(c, []))
    return [" " + norm(x).strip() + " " for x in out if len(norm(x).strip()) >= 3]


def build_profiles(decs, U):
    prof = {}
    for d in decs:
        genes, prots, orgs, missing = set(), set(), set(), []
        for a in d["accessions"]:
            u = U.get(a) or {}
            if not u.get("active"):
                missing.append(f"{a}:{u.get('reason', 'not fetched')}")
                continue
            for g in u["genes"]:
                genes.update(gene_variants(g))
            for p in u["protein_names"]:
                prots.update(protein_variants(p))
                # symbol-like tokens inside protein names (e.g. "eIF-1", "PEN-2") are matched like gene names
                # (not from names of partners: "Hsp90 co-chaperone", "Cdk5 activator"; the letter stem must be a
                # short symbol, not a word such as "factor4" or "oxygenase1")
                if re.search(r"co.?chaperone|binding|interacting|associated|activator|regulator|inhibitor|"
                             r"substrate|partner|anchor|homolog", p, re.I):
                    continue
                for t in norm(p).split():
                    stem = re.match(r"[a-z]*", t).group(0)
                    if len(t) >= 4 and stem and re.search(r"\d", t) and len(stem) <= 6 and \
                            stem not in GENERIC_WORDS and t not in GENERIC_WORDS:
                        genes.update(gene_variants(t))
            orgs.update(organism_variants(u))
        disp_genes = sorted({g for a in d["accessions"] for g in (U.get(a) or {}).get("genes", [])})
        disp_prot = next(((U.get(a) or {}).get("protein_names", [""])[0] for a in d["accessions"]
                          if (U.get(a) or {}).get("active") and (U.get(a) or {}).get("protein_names")), "")
        disp_org = next(((U.get(a) or {}).get("organism") for a in d["accessions"] if (U.get(a) or {}).get("active")), "")
        prof[d["decision_id"]] = {"genes": genes, "prots": prots, "orgs": orgs, "accs": set(d["accessions"]),
                                  "missing": missing, "disp_genes": disp_genes, "disp_prot": disp_prot,
                                  "disp_org": disp_org}
    return prof


GENE_LIKE = re.compile(r"\b(?=[A-Za-z0-9-]*[A-Za-z])(?=[A-Za-z0-9-]*\d)[A-Za-z][A-Za-z0-9-]{1,}\b|\b[A-Z][A-Z]{2,}[a-z]?\b")
BINOMIAL = re.compile(r"\b([A-Z][a-z]{3,}) ([a-z]{4,})\b")


def studied_terms(studied):
    terms = set()
    for m in GENE_LIKE.finditer(studied or ""):
        t = norm(m.group(0)).strip()
        if len(t.replace(" ", "")) >= 3 and t not in GENE_STOP and not re.fullmatch(r"(19|20)\d\d", t) and \
                t not in ("pdb", "gst", "his6", "dna", "rna", "atp", "nmr", "cryo", "em", "fret", "itc"):
            terms.add(" " + t + " ")
    for m in BINOMIAL.finditer(studied or ""):
        terms.add(" " + norm(m.group(0)).strip() + " ")
    return terms


def title_similarity(a, b):
    ta = set(norm(a).split()) - GENERIC_WORDS
    tb = set(norm(b).split()) - GENERIC_WORDS
    if not ta or not tb:
        return None
    return len(ta & tb) / min(len(ta), len(tb))


def claimed_title_core(t):
    """Reviewer 'title' fields are often 'Author et al. (2010), short description' or database labels.
    Keep only text that looks like a real title: strip an author/year prefix, require >= 5 distinctive words."""
    if not t or re.search(r"\b(PMID|PDB|UniProt|RCSB|NCBI|InterPro|PROSITE|doi)\b", t, re.I):
        return None
    t = re.sub(r"^[^,:]{0,60}?(et al\.?|and [A-Z][a-z]+)\s*\(?(19|20)\d\d[a-z]?\)?\s*[,:.]?\s*", "", t)
    t = re.sub(r"\s*\((19|20)\d\d\)\s*$", "", t)
    if len(set(norm(t).split()) - GENERIC_WORDS) < 5:
        return None
    return t


def lit_keys(c, rec):
    """Identifiers a literature record is known by (for checking that IDs given in one source agree)."""
    keys = set()
    if c["id_type"] not in ("pmid", "pmcid", "doi") or rec.get("resolves") != "yes":
        return []
    keys.add(f"{c['id_type']}:{c['id'].lower()}")
    for k in ("pmid", "pmcid", "doi"):
        v = rec.get(k)
        if v:
            keys.add(f"{k}:{str(v).lower()}")
    return sorted(keys)


def match_one(c, rec, prof):
    """-> dict(match_level, matched_on, where)"""
    title, abstract, extra = rec.get("title") or "", rec.get("abstract") or "", rec.get("extra") or ""
    page = rec.get("page_text") or ""
    fields = [("title", title), ("abstract", abstract), ("mesh/keywords", extra), ("page", page)]
    # accession
    if c["id_type"] == "uniprot":
        if c["id"] in prof["accs"]:
            return {"level": "accession", "on": c["id"], "where": "cited UniProt accession is the decision's"}
    if c["id_type"] == "pdb" and set(rec.get("uniprot_ids") or []) & prof["accs"]:
        return {"level": "accession", "on": ",".join(sorted(set(rec["uniprot_ids"]) & prof["accs"])),
                "where": "PDB entity UniProt mapping"}
    for where, t in fields:
        for a in prof["accs"]:
            if t and re.search(r"\b" + re.escape(a) + r"\b", t):
                return {"level": "accession", "on": a, "where": where}
    org_hit = None
    for where, t in fields:
        if not t:
            continue
        nt = norm(t)
        for g in prof["genes"]:
            if g in nt:
                return {"level": "gene_or_protein", "on": g.strip(), "where": where,
                        "organism_also": any(o in norm(title + " " + abstract + " " + extra + " " + page)
                                             for o in prof["orgs"])}
        for p in prof["prots"]:
            if p in nt:
                return {"level": "gene_or_protein", "on": p.strip(), "where": where,
                        "organism_also": any(o in norm(title + " " + abstract + " " + extra + " " + page)
                                             for o in prof["orgs"])}
        if org_hit is None:
            for o in prof["orgs"]:
                if o in nt:
                    org_hit = {"level": "organism_only", "on": o.strip(), "where": where}
                    break
    return org_hit or {"level": "none", "on": "", "where": ""}


def stage_match():
    cites = load_jsonl("citations.jsonl")
    decs = load_jsonl("decisions.jsonl")
    with gzip.open(os.path.join(ARGS.out_tmp, "resolved.json.gz"), "rt", encoding="utf-8") as f:
        R = json.load(f)
    with gzip.open(os.path.join(ARGS.out_tmp, "uniprot_decisions.json.gz"), "rt", encoding="utf-8") as f:
        U = json.load(f)
    decs = [d for d in decs if all(a in U for a in d["accessions"])] if ARGS.limit else decs
    dset = {d["decision_id"] for d in decs}
    cites = [c for c in cites if c["decision_id"] in dset and c["id"] in R.get(c["id_type"], {})]
    prof = build_profiles(decs, U)
    out = []
    for c in cites:
        rec = R[c["id_type"]][c["id"]]
        p = prof[c["decision_id"]]
        row = dict(c)
        row["resolves"] = rec["resolves"]
        row["how"] = rec.get("how", "")
        row["record_title"] = (rec.get("title") or "")[:300]
        has_text = bool((rec.get("title") or "").strip() or (rec.get("abstract") or "").strip() or
                        (rec.get("page_text") or "").strip() or rec.get("uniprot_ids") or c["id_type"] == "uniprot")
        row["has_abstract"] = bool((rec.get("abstract") or "").strip())
        row["text_basis"] = ("page" if rec.get("page_text") else
                             "title+abstract" if row["has_abstract"] else "title only" if rec.get("title") else "none")
        if rec["resolves"] == "yes" and has_text:
            m = match_one(c, rec, p)
        elif rec["resolves"] == "yes":
            m = {"level": "no_text", "on": "", "where": ""}
        else:
            m = {"level": "not_checked", "on": "", "where": ""}
        if c["id_type"] == "uniprot" and rec.get("uniprot") and m["level"] in ("none", "organism_only"):
            # database record: also count gene overlap via UniProt gene names of the cited entry
            cg = {v for g in rec["uniprot"]["genes"] for v in gene_variants(g)}
            if cg & p["genes"]:
                m = {"level": "gene_or_protein", "on": sorted(cg & p["genes"])[0].strip(),
                     "where": "cited UniProt entry gene names", "organism_also": False}
        row["match_level"] = m["level"]
        row["matched_on"] = m["on"]
        row["matched_where"] = m["where"]
        row["organism_also"] = m.get("organism_also")
        st = studied_terms(c.get("studied"))
        alltext = norm(" ".join(str(rec.get(k) or "") for k in ("title", "abstract", "extra", "page_text")))
        row["studied_entity_match"] = (any(t in alltext for t in st) if st and m["level"] not in ("not_checked", "no_text") else None)
        row["title_similarity"] = None
        ct = claimed_title_core(c.get("claimed_title"))
        if ct and c["id_type"] in ("pmid", "pmcid", "doi") and rec.get("title"):
            ts = title_similarity(ct, rec["title"])
            row["title_similarity"] = None if ts is None else round(ts, 2)
        row["lit_keys"] = lit_keys(c, rec)
        ry = rec.get("year")
        row["year_mismatch"] = (abs(int(ry) - c["claimed_year"]) > 1) if (ry and c.get("claimed_year") and
                                                                         c["id_type"] in ("pmid", "pmcid", "doi")) else None
        row["protein"] = p["disp_prot"]
        row["genes"] = ",".join(p["disp_genes"][:8])
        row["organism"] = p["disp_org"]
        row["accessions"] = ",".join(sorted(p["accs"]))
        row["uniprot_missing"] = ";".join(p["missing"])
        out.append(row)
    # identifiers given inside one source (e.g. a PMC URL plus 'PMID ...' in its reference text) must name one paper
    by_src = collections.defaultdict(list)
    for r in out:
        if r["lit_keys"]:
            by_src[(r["decision_id"], r["origin"], r["source_index"])].append(r)
    for rs in by_src.values():
        for r in rs:
            r["same_source_conflict"] = None
        if len(rs) < 2:
            continue
        for r in rs:
            others = [o for o in rs if o is not r]
            ks = set(r["lit_keys"])
            # conflict only if the other record's identifiers are fully known and disjoint
            conflict = [o["id"] for o in others if not (ks & set(o["lit_keys"])) and
                        any(k.startswith("pmid:") for k in ks) and any(k.startswith("pmid:") for k in o["lit_keys"])]
            r["same_source_conflict"] = ",".join(conflict) if conflict else False
    with gzip.open(os.path.join(ARGS.out_tmp, "matched.jsonl.gz"), "wt", encoding="utf-8") as f:
        for r in out:
            f.write(json.dumps(r) + "\n")
    log("match levels:", collections.Counter(r["match_level"] for r in out))


# ----------------------------------------------------------------------------- report
TITLE_FLAG = 0.3
CATS = ["accepted_experimental", "accepted_inference", "rejected"]


def stage_report():
    decs = load_jsonl("decisions.jsonl")
    with gzip.open(os.path.join(ARGS.out_tmp, "matched.jsonl.gz"), "rt", encoding="utf-8") as f:
        rows = [json.loads(l) for l in f]
    os.makedirs(ARGS.out, exist_ok=True)
    S = {"generated": time.strftime("%Y-%m-%d %H:%M %Z"), "inputs": EXPECTED_SHA}

    def frac(a, b):
        return round(a / b, 4) if b else None

    def summarize(rs, ds):
        cited_decs = {r["decision_id"] for r in rs}
        res = collections.Counter(r["resolves"] for r in rs)
        ml = collections.Counter(r["match_level"] for r in rs)
        n = len(rs)
        checkable = [r for r in rs if r["resolves"] == "yes"]
        textful = [r for r in checkable if r["match_level"] not in ("no_text",)]
        prot = sum(r["match_level"] in ("accession", "gene_or_protein") for r in textful)
        org = sum(r["match_level"] == "organism_only" for r in textful)
        # decision level
        bydec = collections.defaultdict(list)
        for r in rs:
            bydec[r["decision_id"]].append(r)
        all_res = sum(all(x["resolves"] == "yes" for x in v) for v in bydec.values())
        any_bad = sum(any(x["resolves"] in ("no",) for x in v) for v in bydec.values())
        any_prot = sum(any(x["match_level"] in ("accession", "gene_or_protein") for x in v) for v in bydec.values())
        any_prot_or_org = sum(any(x["match_level"] in ("accession", "gene_or_protein", "organism_only") for x in v)
                              for v in bydec.values())
        return {"decisions": len(ds), "decisions_citing": len(cited_decs),
                "citation_instances": n, "unique_identifiers": len({(r["id_type"], r["id"]) for r in rs}),
                "by_id_type": dict(collections.Counter(r["id_type"] for r in rs)),
                "resolves": dict(res), "frac_resolves": frac(res["yes"], n),
                "frac_resolves_excl_blocked": frac(res["yes"], n - res["blocked"] - res["unreachable"]),
                "frac_does_not_resolve": frac(res["no"] + res["soft404?"], n),
                "match_levels": dict(ml),
                "resolved_with_text": len(textful),
                "frac_protein_match_of_resolved_with_text": frac(prot, len(textful)),
                "frac_organism_only_of_resolved_with_text": frac(org, len(textful)),
                "frac_no_match_of_resolved_with_text": frac(len(textful) - prot - org, len(textful)),
                "gene_match_without_organism": sum(r["match_level"] == "gene_or_protein" and r.get("organism_also") is False
                                                   for r in textful),
                "decisions_all_citations_resolve": all_res,
                "decisions_with_a_nonresolving_citation": any_bad,
                "decisions_with_protein_matching_citation": any_prot,
                "frac_decisions_with_protein_matching_citation": frac(any_prot, len(bydec)),
                "frac_decisions_with_protein_or_organism_matching_citation": frac(any_prot_or_org, len(bydec)),
                "no_match_but_matches_reviewer_stated_studied_entity":
                    sum(r["match_level"] == "none" and r["studied_entity_match"] is True for r in textful),
                "claimed_title_checked": sum(r["title_similarity"] is not None for r in rs),
                "claimed_title_low_similarity(<0.3)": sum(r["title_similarity"] is not None and r["title_similarity"] < TITLE_FLAG
                                                          for r in rs),
                "same_source_id_conflicts": sum(bool(r.get("same_source_conflict")) for r in rs),
                "year_checked": sum(r["year_mismatch"] is not None for r in rs),
                "year_mismatch(>1y)": sum(bool(r["year_mismatch"]) for r in rs)}
    for split in ("test", "dev", "combined"):
        S[split] = {}
        for cat in CATS + ["all"]:
            ds = [d for d in decs if (split == "combined" or d["split"] == split) and (cat == "all" or d["category"] == cat)]
            dset = {d["decision_id"] for d in ds}
            rs = [r for r in rows if r["decision_id"] in dset]
            S[split][cat] = summarize(rs, ds)
    S["host_class_resolves"] = {hc: dict(collections.Counter(r["resolves"] for r in rows if r["host_class"] == hc))
                                for hc in sorted({r["host_class"] for r in rows})}
    S["host_class_match"] = {hc: dict(collections.Counter(r["match_level"] for r in rows if r["host_class"] == hc))
                             for hc in sorted({r["host_class"] for r in rows})}
    S["blocked_hosts"] = dict(collections.Counter(host_of(r["id"]) for r in rows if r["resolves"] in ("blocked", "unreachable")).most_common(40))
    S["unique_ids_by_resolution"] = dict(collections.Counter(v for v in {(r["id_type"], r["id"]): r["resolves"] for r in rows}.values()))
    S["extract_notes"] = json.load(open(os.path.join(ARGS.out_tmp, "extract_notes.json")))
    S["decisions_with_inactive_uniprot"] = len({r["decision_id"] for r in rows if r["uniprot_missing"]})
    json.dump(S, open(os.path.join(ARGS.out, "summary.json"), "w"), indent=1)

    # failures TSV
    def ftype(r):
        if r["resolves"] in ("no", "soft404?"):
            return "does_not_resolve" if r["resolves"] == "no" else "does_not_resolve(soft-404 page)"
        if r["resolves"] in ("blocked",):
            return "blocked"
        if r["resolves"] == "unreachable":
            return "unreachable(timeout/5xx/connection)"
        if r["resolves"] == "unknown":
            return "resolution_unknown"
        if r.get("same_source_conflict"):
            return "ids_in_same_source_name_different_papers"
        if r["match_level"] == "none":
            return ("resolves_no_match(but_matches_reviewer_stated_studied_entity)" if r["studied_entity_match"]
                    else "resolves_no_protein_gene_organism_match")
        if r["match_level"] == "no_text":
            return "resolves_no_text_to_match"
        if r["title_similarity"] is not None and r["title_similarity"] < TITLE_FLAG:
            return "claimed_title_differs(inspect)"
        return None
    cols = ["decision_id", "split", "category", "id_type", "identifier", "host_class", "failure_type", "how",
            "record_title", "claimed_title", "title_similarity", "studied_entity(reviewer)", "studied_entity_match",
            "match_level", "protein", "genes", "organism", "accessions", "source_section", "raw_url"]
    nfail = collections.Counter()
    with open(os.path.join(ARGS.out, "failures.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["split"], r["category"], r["decision_id"])):
            t = ftype(r)
            if not t:
                continue
            nfail[t] += 1
            vals = [r["decision_id"], r["split"], r["category"], r["id_type"], r["id"], r["host_class"], t, r["how"],
                    r["record_title"], r.get("claimed_title") or "", r["title_similarity"], r.get("studied") or "",
                    r["studied_entity_match"], r["match_level"], r["protein"], r["genes"], r["organism"],
                    r["accessions"], str(r.get("section") or ""), r.get("raw_url") or ""]
            f.write("\t".join(re.sub(r"[\t\r\n]+", " ", "" if v is None else str(v)) for v in vals) + "\n")
    S["failures_by_type"] = dict(nfail)
    json.dump(S, open(os.path.join(ARGS.out, "summary.json"), "w"), indent=1)
    # full per-citation results (compact)
    with gzip.open(os.path.join(ARGS.out, "citation_results.tsv.gz"), "wt", encoding="utf-8") as f:
        cc = ["decision_id", "split", "category", "origin", "id_type", "id", "host_class", "resolves", "match_level",
              "matched_on", "matched_where", "organism_also", "studied_entity_match", "text_basis", "title_similarity",
              "year_mismatch", "record_title"]
        f.write("\t".join(cc) + "\n")
        for r in rows:
            f.write("\t".join(re.sub(r"[\t\r\n]+", " ", "" if r.get(k) is None else str(r.get(k))) for k in cc) + "\n")
    # decisions in which no cited record names the protein/gene (the per-decision view of the failures)
    bydec = collections.defaultdict(list)
    for r in rows:
        bydec[r["decision_id"]].append(r)
    nd = collections.Counter()
    with open(os.path.join(ARGS.out, "decisions_without_protein_match.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("decision_id\tsplit\tcategory\tprotein\tgenes\torganism\tn_citations\tn_resolve\tbest_level\t"
                "any_matches_reviewer_stated_entity\tidentifiers\n")
        order = ["accession", "gene_or_protein", "organism_only", "none", "no_text", "not_checked"]
        for did, rs in sorted(bydec.items()):
            levels = {x["match_level"] for x in rs}
            if levels & {"accession", "gene_or_protein"}:
                continue
            best = next(l for l in order if l in levels)
            st = any(x["studied_entity_match"] for x in rs)
            nd[(rs[0]["split"], rs[0]["category"], best, st)] += 1
            f.write("\t".join(re.sub(r"[\t\r\n]+", " ", str(v)) for v in [
                did, rs[0]["split"], rs[0]["category"], rs[0]["protein"], rs[0]["genes"], rs[0]["organism"], len(rs),
                sum(x["resolves"] == "yes" for x in rs), best, st, " ".join(x["id"] for x in rs)]) + "\n")
    S["decisions_without_protein_match"] = {"|".join(map(str, k)): v for k, v in sorted(nd.items())}
    json.dump(S, open(os.path.join(ARGS.out, "summary.json"), "w"), indent=1)
    log("failures:", dict(nfail))
    write_md(S)


def write_md(S):
    L = []
    names = {"accepted_experimental": "Accepted, experimental", "accepted_inference": "Accepted, inference",
             "rejected": "Rejected", "all": "All"}
    for split in ("test", "dev", "combined"):
        L.append(f"\n### {split}\n")
        L.append("| Category | Decisions | Citing | Citations | Unique IDs | Resolve | Not resolve | Blocked/unreach. | "
                 "Protein/gene match* | Organism only* | No match* | Decisions w/ protein-matching citation |")
        L.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for cat in CATS + ["all"]:
            s = S[split][cat]
            r = s["resolves"]
            pct = lambda x: "" if x is None else f"{100*x:.1f}%"
            L.append(f"| {names[cat]} | {s['decisions']} | {s['decisions_citing']} | {s['citation_instances']} | "
                     f"{s['unique_identifiers']} | {pct(s['frac_resolves'])} | {pct(s['frac_does_not_resolve'])} | "
                     f"{r.get('blocked', 0) + r.get('unreachable', 0)} | {pct(s['frac_protein_match_of_resolved_with_text'])} | "
                     f"{pct(s['frac_organism_only_of_resolved_with_text'])} | {pct(s['frac_no_match_of_resolved_with_text'])} | "
                     f"{pct(s['frac_decisions_with_protein_matching_citation'])} |")
    L.append("\n*Of citations that resolve and have text to match.\n")
    open(os.path.join(ARGS.out, "summary.md"), "w", encoding="utf-8").write("\n".join(L))


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--release", help="the downloaded ArchipelaGO dataset")
    source.add_argument("--data", help="a folder with the review records (test_ and dev_resolved_reviews.json)")
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-tmp", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("stages", nargs="+")
    ARGS = ap.parse_args()
    ARGS.out_tmp = ARGS.out_tmp or ARGS.cache
    os.makedirs(ARGS.cache, exist_ok=True)
    os.makedirs(ARGS.out_tmp, exist_ok=True)
    for s in ARGS.stages:
        log("== stage", s)
        globals()["stage_" + s]()


if __name__ == "__main__":
    main()
