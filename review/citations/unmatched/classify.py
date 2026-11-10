"""Qualitative follow-up on the 147 label-review decisions whose cited records do not name the protein
(best_level = none and any_matches_reviewer_stated_entity = False in the audit's decisions_without_protein_match.tsv), as
reported in Supplementary Note 3: every cited record exists, most are on an ortholog, paralog, family or database entry as the
decision says, and 3 do not fit the protein or the claim (manual_labels.tsv, category d).

Stages (python -I classify.py --cache DIR --audit DIR [--records DIR] stage ...):
  fetch    -> Europe PMC full-text XML for every cited record that is in PMC (cached under <cache>/unmatched_147)
              PMIDs and DOIs are first mapped to a PMCID with a Europe PMC search. Uses urllib, <= 3 requests/s, generic UA.
  gather   -> dossier.jsonl (everything needed to read a decision) and features.tsv
              (one row per decision: automatic features and a provisional category)
  report   -> reads manual_labels.tsv (the hand labels of the paper's reading), writes decisions.tsv and summary_tables.md

The 147 are those of the audit run on the review records (audit_citations.py --data), which, like gather's reviewer terms
and rationales, are not part of the release; an audit run from the release finds a different set of decisions without a
protein match. manual_labels.tsv holds the paper's hand label of each of the 147, with its confidence and reason.

Automatic features per decision
  name terms   UniProt gene names / protein names (audit's variants), plus 'reviewer terms': gene-like tokens and binomials
               that the reviewer itself wrote in `studied_protein_or_organism`, `established` and the rationale.
  where found  title+abstract, full text body (reference list removed; Europe PMC XML only), page text for URL citations.
Provisional category (rules, applied before hand reading; manual_labels.tsv overrides every row):
  c  every cited record is a database/resource page (InterPro, PROSITE, UniProt, GEO, PDB, genome record...)
  a  a UniProt-derived name (gene, protein name, synonym) is found in the full text body
  b? no UniProt-derived name in the text but reviewer terms are, or the rationale speaks of ortholog/homolog/paralog/family
  d? full text available, nothing found, no homology language
  e  no full text for any literature record and nothing found
"""
import argparse, collections, csv, gzip, importlib.util, json, os, re, sys, time
import urllib.error, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.join(HERE, "..", "audit_citations.py")
UA = "Mozilla/5.0 (compatible; citation-check/1.0)"
ARGS = None


def load_audit():
    spec = importlib.util.spec_from_file_location("audit_citations", AUDIT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


A = load_audit()
norm, gene_variants, protein_variants = A.norm, A.gene_variants, A.protein_variants


# ------------------------------------------------------------------ inputs
def paths():
    return {
        "cache": ARGS.cache,
        "ft": os.path.join(ARGS.cache, "unmatched_147"),
        "tsv": os.path.join(ARGS.audit, "decisions_without_protein_match.tsv"),
        "test": os.path.join(ARGS.records or "", "test_resolved_reviews.json"),
        "dev": os.path.join(ARGS.records or "", "dev_resolved_reviews.json"),
    }


def load_base():
    P = paths()
    rows = [r for r in csv.DictReader(open(P["tsv"], encoding="utf8"), delimiter="\t")
            if r["best_level"] == "none" and r["any_matches_reviewer_stated_entity"] == "False"]
    assert len(rows) == 147, len(rows)
    ids = {r["decision_id"] for r in rows}
    c = P["cache"]
    decs = {}
    for l in open(os.path.join(c, "decisions.jsonl"), encoding="utf8"):
        d = json.loads(l)
        if d["decision_id"] in ids:
            decs[d["decision_id"]] = d
    cits = collections.defaultdict(list)
    for l in open(os.path.join(c, "citations.jsonl"), encoding="utf8"):
        x = json.loads(l)
        if x["decision_id"] in ids:
            cits[x["decision_id"]].append(x)
    R = json.load(gzip.open(os.path.join(c, "resolved.json.gz"), "rt", encoding="utf8"))
    U = json.load(gzip.open(os.path.join(c, "uniprot_decisions.json.gz"), "rt", encoding="utf8"))
    rev = {}
    for sp in ("test", "dev"):
        for r in json.load(open(P[sp], encoding="utf8")):
            rev[sp + ":" + r["assertion_id"]] = r
    return rows, decs, cits, R, U, rev


def record_for(R, x):
    t = x["id_type"]
    key = x["id"]
    return R.get(t, {}).get(key) or {}


# ------------------------------------------------------------------ Europe PMC
_last = [0.0]


def get(url, timeout=60):
    gap = 0.34 - (time.time() - _last[0])
    if gap > 0:
        time.sleep(gap)
    _last[0] = time.time()
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code in (404, 403, 410):
                return e.code, b""
            time.sleep(2 * (k + 1))
        except Exception:
            time.sleep(2 * (k + 1))
    return 0, b""


def cached(path, fn):
    if os.path.exists(path):
        return open(path, "rb").read()
    b = fn()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "wb").write(b)
    return b


def to_pmcid(id_type, ident, rec):
    """PMCID for a cited record, or ''."""
    if id_type == "pmcid":
        return ident
    q = None
    if id_type == "pmid":
        q = f"EXT_ID:{ident} AND SRC:MED"
    elif id_type == "doi":
        q = f'DOI:"{ident}"'
    elif rec.get("pmid"):
        q = f"EXT_ID:{rec['pmid']} AND SRC:MED"
    if not q:
        return ""
    P = paths()
    fn = os.path.join(P["ft"], "map_" + re.sub(r"[^A-Za-z0-9]", "_", f"{id_type}_{ident}") + ".json")
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=lite&query="
           + urllib.parse.quote(q))
    b = cached(fn, lambda: get(url)[1])
    try:
        js = json.loads(b)
        for r in js["resultList"]["result"]:
            if r.get("pmcid"):
                return r["pmcid"]
    except Exception:
        pass
    return ""


def fulltext_xml(pmcid):
    P = paths()
    fn = os.path.join(P["ft"], pmcid + ".xml")
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    b = cached(fn, lambda: get(url)[1])
    if b:
        return b
    # not in the Europe PMC open-access set: NCBI PMC efetch (empty body when the publisher bars XML download)
    fn2 = os.path.join(P["ft"], pmcid + ".ncbi.xml")
    url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml"
    return cached(fn2, lambda: get(url2)[1])


def xml_body_text(b):
    """Text of an article without its reference list ('' if not an article)."""
    if not b or b"<body" not in b:
        return ""
    try:
        root = ET.fromstring(b)
    except ET.ParseError:
        return re.sub(r"<[^>]+>", " ", b.decode("utf8", "replace"))
    for par in root.iter():
        for ch in list(par):
            if ch.tag == "ref-list":
                par.remove(ch)
    return " ".join("".join(root.itertext()).split())


def stage_fetch():
    rows, decs, cits, R, U, rev = load_base()
    seen = {}
    for did, xs in cits.items():
        for x in xs:
            if x["id_type"] in ("pmid", "pmcid", "doi", "url"):
                rec = record_for(R, x)
                key = (x["id_type"], x["id"])
                if key in seen:
                    continue
                if x["host_class"] == "other_database":
                    continue
                pm = to_pmcid(x["id_type"], x["id"], rec)
                seen[key] = pm
                if pm:
                    fulltext_xml(pm)
    print("records:", len(seen), "with PMCID:", sum(1 for v in seen.values() if v))
    json.dump({f"{a}|{b}": v for (a, b), v in seen.items()},
              open(os.path.join(paths()["ft"], "pmcid_map.json"), "w"), indent=0)


def interpro_name(url):
    """Name of an InterPro entry (the audit's page fetch got nothing from these JavaScript pages); cached; '' otherwise."""
    m = re.search(r"/entry/InterPro/(IPR\d+)", url)
    if not m:
        return ""
    fn = os.path.join(paths()["ft"], "interpro_" + m.group(1) + ".json")
    b = cached(fn, lambda: get(f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{m.group(1)}")[1])
    try:
        md = json.loads(b)["metadata"]
        return f"{md['name']['name']} ({md.get('type')})"
    except Exception:
        return ""


# ------------------------------------------------------------------ terms
STOP_TOKENS = set("go rna dna atp gtp nad nadh coa ec iba iea isa iso ida imp iep igi iss nas tas rcsb pdb nmr pcr cdna mrna trna rrna "
                  "sds page ph uniprot interpro prosite pfam smart goa quickgo pubmed ncbi ensembl blast clustal pdbe cryo em "
                  "gfp hela hek293 hek293t cos7 ko wt n-terminal c-terminal".split())
ACC = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
TOK = re.compile(r"\b(?=[A-Za-z0-9._-]*[A-Za-z])(?=[A-Za-z0-9._-]*(?:\d|[A-Z][A-Za-z]*[A-Z]))[A-Za-z][A-Za-z0-9._-]{1,}\b")
HOMOLOG = re.compile(r"\b(ortholog\w*|homolog\w*|paralog\w*|famil(?:y|ies)|conserved|conservation|domain|superfamily|"
                     r"sequence (?:identity|similarity)|identical|transfer)\b", re.I)


def reviewer_terms(pr):
    texts = []
    for s in pr.get("sources") or []:
        if isinstance(s, dict):
            for k in ("studied_protein_or_organism", "studied_system", "studied_entity", "established"):
                if s.get(k):
                    texts.append(str(s[k]))
    texts.append(pr.get("rationale") or "")
    out = collections.OrderedDict()
    for t in texts:
        for m in TOK.findall(t):
            m = m.strip(".-_")
            n = norm(m).strip()
            if len(n) < 3 or n in STOP_TOKENS or ACC.match(m) or re.match(r"^(go|ipr|pf|ps|pdb|pmid|pmc|doi)\W*\d", m, re.I):
                continue
            if n.isdigit() or re.fullmatch(r"\d+(\.\d+)?[a-z]{0,2}", n):
                continue
            if re.fullmatch(r"(n|c|ph|mg|nm|kda|ec)\d*", n) or n in A.GENE_STOP:
                continue
            out[m] = 1
    return list(out)[:25]


def find_terms(terms, ntext):
    hits = []
    for label, vs in terms:
        for v in vs:
            if v in ntext:
                hits.append(label)
                break
    return hits


def build_terms(d, U, pr):
    """UniProt-derived terms [(label, [variants])] and reviewer terms."""
    uni = []
    for a in d["accessions"]:
        u = U.get(a) or {}
        if not u.get("active"):
            continue
        for g in u["genes"]:
            v = gene_variants(g)
            if v:
                uni.append((g, v))
        for p in u["protein_names"]:
            v = protein_variants(p)
            if v:
                uni.append((p, v))
    seen, uni2 = set(), []
    for l, v in uni:
        if l not in seen:
            seen.add(l)
            uni2.append((l, v))
    rv = []
    for t in reviewer_terms(pr):
        v = [" " + norm(t).strip() + " "]
        rv.append((t, v))
    return uni2, rv


# ------------------------------------------------------------------ gather
def stage_gather():
    rows, decs, cits, R, U, rev = load_base()
    P = paths()
    pmap = {}
    mp = os.path.join(P["ft"], "pmcid_map.json")
    if os.path.exists(mp):
        pmap = {tuple(k.split("|", 1)): v for k, v in json.load(open(mp)).items()}
    dossier, feats = [], []
    for r in rows:
        did = r["decision_id"]
        d = decs[did]
        rv = rev[did]
        pr = rv.get("pro_review") or {}
        gs = rv.get("original_goa_support") or []
        ev = sorted({e for g in gs for e in g.get("evidence_codes", [])})
        uni, rvt = build_terms(d, U, pr)
        srcs = []
        for s in pr.get("sources") or []:
            if isinstance(s, dict):
                srcs.append({"url": s.get("url"), "locator": s.get("section") or s.get("locator") or s.get("location"),
                             "established": s.get("established"),
                             "studied": s.get("studied_protein_or_organism") or s.get("studied_system"),
                             "access": s.get("access")})
        cl = []
        any_ft = False
        lit_hit_uni = lit_hit_rv = False
        n_lit = n_db = 0
        for x in cits[did]:
            rec = record_for(R, x)
            db = x["id_type"] in ("uniprot", "pdb") or x["host_class"] == "other_database"
            dbname = interpro_name(x.get("raw_url") or "") if db else ""
            title, abs_ = rec.get("title") or "", rec.get("abstract") or ""
            ta = norm(title + " " + abs_ + " " + (rec.get("extra") or ""))
            pm = pmap.get((x["id_type"], x["id"]), "")
            ftext = ""
            if pm:
                ftext = xml_body_text(fulltext_xml(pm))
            nft = norm(ftext) if ftext else ""
            hu_ta = find_terms(uni, ta)
            hr_ta = find_terms(rvt, ta)
            hu_ft = find_terms(uni, nft) if nft else []
            hr_ft = find_terms(rvt, nft) if nft else []
            if db:
                n_db += 1
            else:
                n_lit += 1
                if ftext:
                    any_ft = True
                lit_hit_uni |= bool(hu_ta or hu_ft)
                lit_hit_rv |= bool(hr_ta or hr_ft)
            cl.append({"id_type": x["id_type"], "id": x["id"], "raw_url": x.get("raw_url"), "host_class": x["host_class"],
                       "resolves": rec.get("resolves"), "db": db, "db_entry_name": dbname, "title": title, "abstract": abs_[:1500],
                       "pmcid": pm, "fulltext_chars": len(ftext), "claimed_year": x.get("claimed_year"),
                       "hits_uniprot_ta": hu_ta[:8], "hits_reviewer_ta": hr_ta[:8],
                       "hits_uniprot_ft": hu_ft[:8], "hits_reviewer_ft": hr_ft[:8]})
        homol = bool(HOMOLOG.search((pr.get("rationale") or "")))
        if n_lit == 0:
            prov = "c"
        elif lit_hit_uni:
            prov = "a"
        elif lit_hit_rv or homol:
            prov = "b?"
        elif any_ft:
            prov = "d?"
        else:
            prov = "e"
        u0 = next((U[a] for a in d["accessions"] if (U.get(a) or {}).get("active")), {})
        dos = {"decision_id": did, "split": d["split"], "category": d["category"], "go_id": d["go_id"], "go_name": d["go_name"],
               "protein": r["protein"], "genes": r["genes"], "organism": r["organism"], "accessions": d["accessions"],
               "uniprot_names": sorted({p for a in d["accessions"] for p in (U.get(a) or {}).get("protein_names", [])})[:6],
               "evidence_codes": ev, "review_support_type": rv.get("review_support_type"),
               "rationale": pr.get("rationale"), "limitations": pr.get("limitations"), "sources": srcs,
               "citations": cl, "reviewer_terms": [t for t, _ in rvt], "homology_language": homol}
        dossier.append(dos)
        feats.append({"decision_id": did, "split": d["split"], "category": d["category"], "go_id": d["go_id"],
                      "go_name": d["go_name"], "protein": r["protein"], "organism": r["organism"],
                      "evidence_codes": ",".join(ev), "n_literature": n_lit, "n_database": n_db,
                      "fulltext_available": any_ft, "uniprot_name_found": lit_hit_uni,
                      "reviewer_term_found": lit_hit_rv, "homology_language": homol, "provisional": prov,
                      "ids": ";".join(f"{c['id_type']}:{c['id'] if c['id_type'] != 'url' else c['raw_url'][:80]}" for c in cl)})
    with open(os.path.join(HERE, "dossier.jsonl"), "w", encoding="utf8") as f:
        for x in dossier:
            f.write(json.dumps(x, ensure_ascii=False) + "\n")
    with open(os.path.join(HERE, "features.tsv"), "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(feats[0]), delimiter="\t")
        w.writeheader()
        w.writerows(feats)
    print(collections.Counter(x["provisional"] for x in feats))


# ------------------------------------------------------------------ report
def stage_report():
    feats = {r["decision_id"]: r for r in csv.DictReader(open(os.path.join(HERE, "features.tsv"), encoding="utf8"), delimiter="\t")}
    man = {r["decision_id"]: r for r in csv.DictReader(open(os.path.join(HERE, "manual_labels.tsv"), encoding="utf8"), delimiter="\t")}
    assert set(man) == set(feats), (len(set(feats) - set(man)), len(set(man) - set(feats)))
    out = []
    for did, f in feats.items():
        m = man[did]
        out.append({**{k: f[k] for k in ("decision_id", "split", "category", "go_id", "go_name", "protein", "organism",
                                         "evidence_codes", "provisional")},
                    "final_category": m["category"], "confidence": m["confidence"], "reason": m["reason"], "hand_read": m["hand_read"], "ids": f["ids"]})
    with open(os.path.join(HERE, "decisions.tsv"), "w", encoding="utf8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    cats = ["a", "b", "c", "d", "e"]
    lines = []
    def table(title, keyf):
        keys = sorted({keyf(o) for o in out})
        lines.append(f"\n### {title}\n")
        lines.append("| | " + " | ".join(cats) + " | total |")
        lines.append("|---|" + "---:|" * (len(cats) + 1))
        for k in keys:
            sub = [o for o in out if keyf(o) == k]
            lines.append(f"| {k} | " + " | ".join(str(sum(1 for o in sub if o['final_category'] == c)) for c in cats) + f" | {len(sub)} |")
        lines.append("| all | " + " | ".join(str(sum(1 for o in out if o['final_category'] == c)) for c in cats) + f" | {len(out)} |")
    table("By split and decision category", lambda o: f"{o['split']} {o['category']}")
    table("By split", lambda o: o["split"])
    table("By decision category", lambda o: o["category"])
    prov = collections.Counter((o["provisional"], o["final_category"]) for o in out)
    lines.append("\n### Provisional (rule) category vs final (hand) category\n")
    lines.append("| provisional | " + " | ".join(cats) + " |\n|---|" + "---:|" * len(cats))
    for p in sorted({o["provisional"] for o in out}):
        lines.append(f"| {p} | " + " | ".join(str(prov[(p, c)]) for c in cats) + " |")
    open(os.path.join(HERE, "summary_tables.md"), "w", encoding="utf8").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, help="the audit's cache folder (its --out-tmp: decisions.jsonl, citations.jsonl, "
                    "resolved.json.gz, uniprot_decisions.json.gz)")
    ap.add_argument("--audit", required=True, help="the audit's --out folder (decisions_without_protein_match.tsv)")
    ap.add_argument("--records", help="the review records (test_ and dev_resolved_reviews.json; not part of the release)")
    ap.add_argument("stages", nargs="+")
    ARGS = ap.parse_args()
    for s in ARGS.stages:
        globals()["stage_" + s]()


if __name__ == "__main__":
    main()
