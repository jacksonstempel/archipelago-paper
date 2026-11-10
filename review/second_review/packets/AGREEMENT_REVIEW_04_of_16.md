# Protein function review — batch 4 of 16

25 proteins; 25 assertions. Name your output `AGREEMENT_RESULTS_04.json`.

Review the supplied protein–GO molecular-function assertions. Research the evidence for each assertion and decide whether it is sufficiently supported to accept.

Work protein by protein, reusing relevant research across its assertions. The packet contains sequences, identifiers, GO definitions, and existing annotation provenance. Use primary literature and authoritative databases as appropriate. Investigate the underlying evidence rather than treating an existing annotation as confirmation.

## Decision rubric

Return **accept**, **do not accept**, or, exceptionally, **needs adjudication**.

Accept when either:

- **Experimental evidence** establishes the requested function for the protein. Biochemical, structural, interaction, and genetic experiments can qualify; a purified-protein assay is not required. The evidence must establish the particular molecular function, rather than merely a related phenotype or biological process.
- **A well-supported inference** establishes the function. Explain why the experimental evidence from another protein, or a curated functional rule, applies to the supplied protein. Relevant domain conservation, essential residues, orthology, and known functional differences should inform this judgment. Generic family resemblance or a repeated annotation alone is insufficient.

Apply these principles consistently:

- **Match the exact GO definition.** Evidence for a broader or related activity does not establish a more specific claim.
- **Weak activity is not automatically invalid.** Accept convincing, specific activity even if its physiological role is unknown. Do not accept a marginal signal that remains indistinguishable from nonspecific association or experimental artifact.
- **Complex-dependent activity can qualify.** Accept an established functional contribution to a complex; mere complex membership does not suffice. Record whether the evidence supports `enables` or `contributes_to`.
- **Relevant constructs and mature products can qualify.** Experiments on a functional domain or naturally processed product may establish a function of the supplied sequence. Strain, isoform, or construct differences should prevent acceptance only when they materially weaken the attribution.
- **Check the inference, not every historical pathway.** One adequately substantiated route to support is sufficient; you do not need to reconstruct every annotation’s provenance.
- **Failure to find support is not evidence of falsehood.** Return “do not accept” when sufficient support cannot be established. Reserve claims of contradiction for evidence genuinely opposing the assertion.

Be willing to accept justified functional inferences without demanding a new direct assay for every protein. Conversely, before rejecting an assertion, consider whether relevant evidence appears under another identifier, in a functional domain, or in a different study from the one originally cited.

Use **needs adjudication** only when a concrete conflict or ambiguity prevents a defensible decision under this rubric—for example, incompatible authoritative term definitions or strong conflicting evidence. Ordinary missing evidence should result in **do not accept**, not adjudication.

## Deliverable

Return a JSON file containing one entry per supplied assertion:

- `assertion_id`: preserve exactly.
- `decision`: `accept`, `do_not_accept`, or `needs_adjudication`.
- `support_type`: `experimental`, `inference`, or `none_established`.
- `supported_relation`: `enables`, `contributes_to`, or null.
- `rationale`: concise explanation connecting the evidence to the exact GO term.
- `sources`: URLs, relevant figure/section/record, what was established, and the studied protein or organism where relevant. State when only an abstract or limited material was accessible.
- `limitations`: any qualification material to the decision.
- `adjudication_question`: null unless adjudication is requested; explain precisely what remains undecidable.

Distinguish facts you verified from assumptions or inferences. Do not claim sequence comparisons, source access, or measurements you did not perform. Note evidence published after the packet’s annotation snapshot.

Review every supplied assertion. If you cannot finish, list unreviewed assertion IDs separately rather than assigning them an evidence verdict. Include a short summary of completion and decision counts.


## Protein and assertion data

```json
{
 "batch": 4,
 "source_release": "GOA 2026-07-29",
 "proteins": [
  {
   "sequence_sha256": "25ff0a17247e29815e03098be80bf0042c4e2fdae93174e3ee8b98fc321fcffb",
   "sequence": "MALISCPRAPLWLGRLSRRLCARHRACTGTMARPRRFTVELPDCSLTHFVLGDAAGHRDARLTELLGPPGRSYALCVPLVPGEGCGPRVQAARVYQRLQLQLRQGPLQRCELSKLLGYGPGDQASEAQHGFLLRDPSDHPDTRQDLLQLLGSCQEATRLHLAEFQADSQGLLWQRLWELQGDRQVQVDCACVLPAPEPHLHPLLPDLLNSAVFQDRDAARAVLEECTSFIPEAREVLDLVDQCPREIQKGKFQVIAIEGLDATGKTTLTQSVSESLQAVLLQSPPPCISQWRKLFDDEPTIIRRAFYSLGNYLVASEIAKQSAKFPVIVDRYWHSTATYAIATEVSGGLQYLPPAHHPVYQWPRDLLKPDLVLLLTVNSEERVRRLQGRGQEKTKEEAELEVNSVFRQKVEATYQRMENPSCRLVDASPSRETVLQTVLQLIQRAGH",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "D3ZC63"
   ],
   "assertions": [
    {
     "assertion_id": "25ff0a17247e29815e03098be80bf0042c4e2fdae93174e3ee8b98fc321fcffb:GO:0033862",
     "go_id": "GO:0033862",
     "go_name": "UMP kinase activity",
     "go_definition": "\"Catalysis of the reaction: ATP + UMP = ADP + UDP.\" [EC:2.7.4.22]",
     "supporting_assertions": [
      {
       "accession": "D3ZC63",
       "symbol": "Cmpk2",
       "relation": "enables",
       "go_id": "GO:0033862",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q5EBM0"
       ],
       "object_name": "UMP-CMP kinase 2, mitochondrial",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20100610",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/Q5EBM0/entry"
       ]
      },
      {
       "accession": "D3ZC63",
       "symbol": "Cmpk2",
       "relation": "enables",
       "go_id": "GO:0033862",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q5EBM0",
        "ensembl:ENSP00000256722"
       ],
       "object_name": "UMP-CMP kinase 2, mitochondrial",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q5EBM0/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "3f6a3bac413527b7bac6f2a6983780f653985da2f4d12fc7aed047acc0d7b6a9",
   "sequence": "MRNSLIHHTRNLSATKQITGIVCDHVDVFTWGHSPPEQHGNKFSSGKATAPFYNAGQRRCTLCDDRLETSFSAHCGYVGHVARVGILERAMEILQKGEKRAANGKSGGLGGKITDLDESIKALVGTWWSRLNDTKREPALDYKRIPSLSASTTKKRLWRMRFLLQYLRDRGITRYSLTPAKVAGAGDNAFVRSARFERSEMIGDNIVKVVVPDRLVRLFPADEGGVTYKLASIQQLLDSNEGLLEIYDYLGLNNIIGVRLPNNKTKSDVVEALFGELQTFLWASELSCGSCQYPAFPTAEHRYVRALVDHLLNELTHMVIMWRVGSTLENSKEYLAKYLLQGARQNGSFSTSASGATIVKEVDCDRSRYAVLPLLLTFPYSTSSANKSAGGKQPPLQSSSLAAAKSVQERGTPRMLFTSPVAMPAALHSHYNPPLRIVRTIKDYQREILTVLEEEEVVSALIPCTDSSWGMKEGYMNSRETGRRKPPHGQIRAPWHELALEKYPRWRCDTEEAVVRQALELRLMQRKVYTKGTNWRALTTAMAQRLLPNSPLEVDPVGVRPVRQTGRLLLTFEGTTERRAALATSMCFPSLAQSAG",
   "taxon_ids": [
    "5691"
   ],
   "accessions": [
    "Q86MW2"
   ],
   "assertions": [
    {
     "assertion_id": "3f6a3bac413527b7bac6f2a6983780f653985da2f4d12fc7aed047acc0d7b6a9:GO:0004525",
     "go_id": "GO:0004525",
     "go_name": "ribonuclease III activity",
     "go_definition": "\"Catalysis of the endonucleolytic cleavage of RNA with 5'-phosphomonoesters and 3'-OH termini; makes two staggered cuts in both strands of dsRNA, leaving a 3' overhang of 2 nt.\" [PMID:11157775, PMID:15242644]",
     "supporting_assertions": [
      {
       "accession": "Q86MW2",
       "symbol": "MP67",
       "relation": "enables",
       "go_id": "GO:0004525",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR036389"
       ],
       "object_name": "RNA editing complex protein MP67",
       "taxon": [
        "taxon:5691"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036389/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "46c3d2c34b5424b06c500f77429683e454e99639b8763c4b6365844947ca4cae",
   "sequence": "MKAQGETEDSERLSKMSSLLERLHAKFNQNRPWSETIKLVRQVMEKRVVMSSGGHQHLVSCLETLQKALKVTSLPAMTDRLESIARQNGLGSHLSASGTECYITSDMFYVEVQLDPAGQLCDVKVAHHGENPVSCPELVQQLREKNFEEFSKHLKGLVNLYNLPGDNKLKTKMYLALQSLEQDLSKMAIMYWKATNAAPLDKILHGSVGYLTPRSGGHLMNMKYYASPSDLLDDKTASPIILHEKNVPRSLGMNASVTIEGTSAMYKLPIAPLIMGSHPADNKWTPSFSAVTSANSVDLPACFFLKFPQPIPVSKAFVQKLQNCTGIPLFETPPTYLPLYELITQFELSKDPDPLPLNHNMRFYAALPGQQHCYFLNKDAPLPDGQSLQGTLVSKITFQHPGRVPLILNMIRHQVAYNTLIGSCVKRTILKEDSPGLLQFEVCPLSESRFSVSFQHPVNDSLVCVVMDVQDSTHVSCKLYKGLSDALICTDDFIAKVVQRCMSIPVTMRAIRRKAETIQADTPALSLIAETVEDMVKKNLPPASSPGYGMTTGNNPMSGTTTPTNTFPGGPITTLFNMSMSIKDRHESVGHGEDFSKVSQNPILTSLLQITGNGGSTIGSSPTPPHHTPPPVSSMAGNTKNHPMLMNLLKDNPAQDFSTLYGSSPLERQNSSSGSPRMEMCSGSNKAKKKKSSRVPPDKPKHQTEDDFQRELFSMDVDSQNPMFDVSMTADALDTPHITPAPSQCSTPPATYPQPVSHPQPSIQRMVRLSSSDSIGPDVTDILSDIAEEASKLPSTSDDCPPIGTPVRDSSSSGHSQSALFDSDVFQTNNNENPYTDPADLIADAAGSPNSDSPTNHFFPDGVDFNPDLLNSQSQSGFGEEYFDESSQSGDNDDFKGFASQALNTLGMPMLGGDNGEPKFKGSSQADTVDFSIISVAGKALGAADLMEHHSGSQSPLLTTGELGKEKTQKRVKEGNGTGASSGSGPGSDSKPGKRSRTPSNDGKSKDKPPKRKKADTEGKSPSHSSSNRPFTPPTSTGGSKSPGSSGRSQTPPGVATPPIPKITIQIPKGTVMVGKPSSHSQYTSSGSVSSSGSKSHHSHSSSSSSLASASTSGKVKSSKSEGSSSSKLSGSMYASQGSSGSSQSKNSSQTGGKPGSSPITKHGLSSGSSSTKMKPQGKPSSLMNPSISKPNISPSHSRPPGGSDKLASPMKPVPGTPPSSKAKSPISSGSSGSHVSGTSSSSGMKSSSGSASSGSVSQKTPPASNSCTPSSSSFSSSGSSMSSSQNQHGSSKGKSPSRNKKPSLTAVIDKLKHGVVTSGPGGEDPIDSQMGASTNSSNHPMSSKHNTSGGEFQSKREKSDKDKSKVSASGGSVDSSKKTSESKNVGSTGVAKIIISKHDGGSPSIKAKVTLQKPGESGGDGLRPQIASSKNYGSPLISGSTPKHERGSPSHSKSPAYTPQNVDSESESGSSIAERSYQNSPSSEDGIRPLPEYSTEKHKKHKKEKKKVRDKDRDKKKSHSMKPENWSKSPISSDPTASVTNNPILSADRPSRLSPDFMIGEEDDDLMDVALIGN",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q925J9"
   ],
   "assertions": [
    {
     "assertion_id": "46c3d2c34b5424b06c500f77429683e454e99639b8763c4b6365844947ca4cae:GO:0000978",
     "go_id": "GO:0000978",
     "go_name": "RNA polymerase II cis-regulatory region sequence-specific DNA binding",
     "go_definition": "\"Binding to a specific upstream regulatory DNA sequence (transcription factor recognition sequence or binding site) located in cis relative to the transcription start site (i.e., on the same strand of DNA) of a gene transcribed by RNA polymerase II.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "Q925J9",
       "symbol": "Med1",
       "relation": "enables",
       "go_id": "GO:0000978",
       "references": [
        "GO_REF:0000096"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "RGD:1559552"
       ],
       "object_name": "Mediator of RNA polymerase II transcription subunit 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000096.html"
       ]
      },
      {
       "accession": "Q925J9",
       "symbol": "Med1",
       "relation": "enables",
       "go_id": "GO:0000978",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:D3ZRN2",
        "ensembl:ENSRNOP00000059848"
       ],
       "object_name": "Mediator of RNA polymerase II transcription subunit 1",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/D3ZRN2/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "11bf09525931d3daf28901e9c47d7e9f8ca9e8e6864f5a753d6e0fd037b2b52b",
   "sequence": "MKRPSFLPVLIGTGFGSGFSPFAPGTAGALLASIIWIALYFLLPFTALLWTTAALVVLFTFAGIWAANKLESCWGEDPSRVVVDEMVGVWIPLLAVPDNDRWYWYVIAAFALFRIFDIVKPLGVRKMENFKGGVGVMMDDVLAGVYSFILIAVARWVIG",
   "taxon_ids": [
    "226186"
   ],
   "accessions": [
    "Q8A7J9"
   ],
   "assertions": [
    {
     "assertion_id": "11bf09525931d3daf28901e9c47d7e9f8ca9e8e6864f5a753d6e0fd037b2b52b:GO:0008962",
     "go_id": "GO:0008962",
     "go_name": "phosphatidylglycerophosphatase activity",
     "go_definition": "\"Catalysis of the reaction: phosphatidylglycerophosphate + H2O = phosphatidylglycerol + phosphate.\" [PMID:4292860]",
     "supporting_assertions": [
      {
       "accession": "Q8A7J9",
       "symbol": "BT_1525",
       "relation": "enables",
       "go_id": "GO:0008962",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR007686",
        "InterPro:IPR036681"
       ],
       "object_name": "Putative phosphatidylinositol-3-phosphatase",
       "taxon": [
        "taxon:226186"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR007686/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036681/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "9b8972062b132bfb3e39b3ce96daa063bed09dbc4adf406d42e1d7eb87d57331",
   "sequence": "MLRFDFPGVSIGAAHYEEGPTGATVIHIPAGARTAVDARGGAVGLSGGYDFNHAICLAGGAGYGLEAGAGVSGALLERLEYRTGFAELQLVSSAVIYDFSARSTAVYPDKALGRAALEFAVPGEFPQGRAGAGMSASAGKVDWDRTEITGQGAAFRRLGDVRILAVVVPNPVGVIVDRAGTVVRGNYDAQTGVRRHPVFDYQEAFAEQVPPVTEAGNTTISAIVTNVRMSPVELNQFAKQVHSSMHRGIQPFHTDMDGDTLFAVTTDEIDLPTTPGSSRGRLSVNATALGAIASEVMWDAVLEAGK",
   "taxon_ids": [],
   "accessions": [
    "Q14T92"
   ],
   "assertions": [
    {
     "assertion_id": "9b8972062b132bfb3e39b3ce96daa063bed09dbc4adf406d42e1d7eb87d57331:GO:0004177",
     "go_id": "GO:0004177",
     "go_name": "aminopeptidase activity",
     "go_definition": "\"Catalysis of the hydrolysis of a single N-terminal amino acid residue from a polypeptide chain.\" [EC:3.4.11.-, https://www.ebi.ac.uk/merops/about/glossary.shtml#AMINOPEPTIDASE, PMID:24157837]",
     "supporting_assertions": [
      {
       "accession": "Q14T92",
       "symbol": "nylC'",
       "relation": "enables",
       "go_id": "GO:0004177",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN002194607"
       ],
       "object_name": "Putative NylC-analogous protein",
       "taxon": [
        "taxon:388924"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "a9edf59dafc421d7ab3ccedd05eee1d332c1758684b2a215c3573b20be545160",
   "sequence": "MQGRLEDGVVHLPGDARQRFHDSRGYGRPTGGDDLEVAPVEAAHLLSRDDIDGVDGMGLRELLARTGTTLDFVVYKDLRDRGFYLSPAREGWPGVADAADADFLVYPRGKGPWDGEVEHRVRVVGERESIPVSSLGEVVLAIVDEDGDLTYFDTEGDDPEGTAAEDLPADLDAELLSDRALVWDGVDRLYQRGFFGQRLYGRNADSGPLQLSLLEAAYLARADALAIDEADVVSRGRDVEGDRFDRRLAVYAALREAKTVPKSGFKFGSDFRVYTEFESVDDLSHSEFLVRVVAPDHTFVPRDLSLDVRLAGGVRKRMVFALTDDNGEIDWLSVSRLTP",
   "taxon_ids": [
    "309800"
   ],
   "accessions": [
    "O07118"
   ],
   "assertions": [
    {
     "assertion_id": "a9edf59dafc421d7ab3ccedd05eee1d332c1758684b2a215c3573b20be545160:GO:0003676",
     "go_id": "GO:0003676",
     "go_name": "nucleic acid binding",
     "go_definition": "\"Binding to a nucleic acid.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "O07118",
       "symbol": "endA",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR011856"
       ],
       "object_name": "tRNA-splicing endonuclease",
       "taxon": [
        "taxon:309800"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR011856/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "da4c763421c86fa315a07859eb02dd2f66a2e1410595dc8efe716bd3ab70a9f3",
   "sequence": "MKRMEGKLFMVSTKLQAVTKAVLLSTVLSISLLNNEVIKAEQLNMNSQNKYTNFENLKITDKVEDFKEDKEKAKEWGKEKEKEWKLTATEKGKMNNFLDNKNDIKTNYKEITFSMAGSFEDEIKDLKEIDKIFDKANLSSPIITYKNVEPATIGFNKSLTEGNTINSDAMAQFKEQFLDRDIKFDSYLDTHLTVQQVSSKERVILKVKVPSGKGSTTPTKAGIILNNSEYKMLIDNGYMVHVDKVSKVVKKGVECLQVEGTLKKSLDFKNDINAGAHSWGMKNYEEWAKDLTDLQREALDGYARQDYKEINNYLRNQGGNGNEKLDAQIKNISDALGKKPIPENITVYRWCGMPEFGYQISDPLPSLKDFEEQFLNTIKEDKGYMSTSLSSERLAAFGSRKIILRLQVPKGSTGAYLSAIGGFANEKEILLDKDSKYHIDKVTEVIIKGVKRYVVDATLLTN",
   "taxon_ids": [],
   "accessions": [
    "Q844J9"
   ],
   "assertions": [
    {
     "assertion_id": "da4c763421c86fa315a07859eb02dd2f66a2e1410595dc8efe716bd3ab70a9f3:GO:0106274",
     "go_id": "GO:0106274",
     "go_name": "NAD+-protein-arginine ADP-ribosyltransferase activity",
     "go_definition": "\"Catalysis of the reaction: L-arginyl-[protein] + NAD+ = H+ + (ADP-D-ribosyl)-L-arginyl-[protein] + nicotinamide.\" [RHEA:19149]",
     "supporting_assertions": [
      {
       "accession": "Q844J9",
       "symbol": "vip2Ac",
       "relation": "enables",
       "go_id": "GO:0106274",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR001683352"
       ],
       "object_name": "ADP ribosyltransferase domain-containing protein",
       "taxon": [
        "taxon:1428"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "9873c0ca6c6b1a7f08a46f7c0bca7cd3f0a7456fce2c4c8c6cad9067f26dede5",
   "sequence": "MAGTPQPRALGPDALDVSTDDLAGLLAGNTGRIKTVITDQKVIAGIGNAYSDEILHVAKISPFATAGKLSGAQLTCLHEAMASVLSDAVRRSVGQGAAMLKGEKRSGLRVHARTGLPCPVCGDTVREVSFADKSFQYCPTCQTGGKALADRRMSRLLK",
   "taxon_ids": [
    "83332"
   ],
   "accessions": [
    "A0A1R3XWW9",
    "A0ABX2VKC1",
    "A0ACD3TT69",
    "A5U0Y4",
    "L0T864"
   ],
   "assertions": [
    {
     "assertion_id": "9873c0ca6c6b1a7f08a46f7c0bca7cd3f0a7456fce2c4c8c6cad9067f26dede5:GO:0003676",
     "go_id": "GO:0003676",
     "go_name": "nucleic acid binding",
     "go_definition": "\"Binding to a nucleic acid.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "L0T864",
       "symbol": "fpg2",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR010979"
       ],
       "object_name": "Uncharacterized formamidopyrimidine-DNA glycosylase-like protein",
       "taxon": [
        "taxon:83332"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010979/"
       ]
      },
      {
       "accession": "A5U0Y4",
       "symbol": "MRA_0951",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR010979"
       ],
       "object_name": "FPG-type domain-containing protein",
       "taxon": [
        "taxon:419947"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010979/"
       ]
      },
      {
       "accession": "A0ACD3TT69",
       "symbol": "FHI80_06390",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR010979"
       ],
       "object_name": "FPG-type domain-containing protein",
       "taxon": [
        "taxon:1305738"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010979/"
       ]
      },
      {
       "accession": "A0ABX2VKC1",
       "symbol": "A7J32_07705",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR010979"
       ],
       "object_name": "FPG-type domain-containing protein",
       "taxon": [
        "taxon:1844474"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010979/"
       ]
      },
      {
       "accession": "A0A1R3XWW9",
       "symbol": "BQ2027_MB0969",
       "relation": "enables",
       "go_id": "GO:0003676",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR010979"
       ],
       "object_name": "FPG-type domain-containing protein",
       "taxon": [
        "taxon:233413"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR010979/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "1707ef48ea0463788f2996eeb4224269ffe13c8647087be050218ea3e5a16d47",
   "sequence": "MQSPAFAYREHAEHSPRAMYMRDWRSALRTPTYRIGNRTVRFNYHFALLIVCAFVLVLLFYFARQGSHSSSDGYWLRRSYTDARSLDNGGPVVQAYNATYPLTSPMVLRGGVINYRIAMIADLDTSSKVSKGDGSSTWRSYLKKGYLTYTVARSEIQISWDDGAPIVLESAFALKGRGMELSELVTFNGRLLTFDDRTGLIYEIVNDKPIPWVILLDGDGHSAKGFKAEWATVKEQTLYVGSMGKEWTTSAGDFENNNPMYVKAITPSGEVRSLNWVDNFKQLRLQSMQITWPGYMIHESGTWSEERNRWFFLPRRCSKEKYNETKDEHMGCNVLVSADESFTNVETVRLDPENTTPTHGFSSFKFLPGTDDSIIVALKSEELNGKTATFITAFDIAGKTLLPETRIETDYKYEGFEFI",
   "taxon_ids": [
    "7227"
   ],
   "accessions": [
    "Q7K3F9",
    "Q9VGN8"
   ],
   "assertions": [
    {
     "assertion_id": "1707ef48ea0463788f2996eeb4224269ffe13c8647087be050218ea3e5a16d47:GO:0004382",
     "go_id": "GO:0004382",
     "go_name": "GDP phosphatase activity",
     "go_definition": "\"Catalysis of the reaction: GDP + H2O = GMP + phosphate.\" [PMID:2989286, RHEA:22156]",
     "supporting_assertions": [
      {
       "accession": "Q9VGN8",
       "symbol": "Dmel\\CG5276",
       "relation": "enables",
       "go_id": "GO:0004382",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "PANTHER:PTN000989481",
        "UniProtKB:Q8WVQ1",
        "WB:WBGene00017244"
       ],
       "object_name": "Apyrase",
       "taxon": [
        "taxon:7227"
       ],
       "annotation_date": "20170228",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/Q8WVQ1/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "97f334e055048d81c04ba0a6a113c16782448fe9ff0db2dc54f20cfd528f2794",
   "sequence": "MTLAAMTAAASQLDNAAPDDVEITDSSAAAEYIADGCLVDGPLGRVGLEMEAHCFDPADPFRRPSWEEITEVLEWLSPLPGGSVVSVEPGGAVELSGPPADGVLAAIGAMTRDQAVLRSALANAGLGLVFLGADPLRSPVRVNPGARYRAMEQFFAASHSGVPGAAMMTSTAAIQVNLDAGPQEGWAERVRLAHALGPTMIAIAANSPMLGGRFSGWQSTRQRVWGQMDSARCGPILGASGDHPGIDWAKYALKAPVMMVRSPDTQDTRAVTDYVPFTDWVDGRVLLDGRRATVADLVYHLTTLFPPVRPRQWLEIRYLDSVPDEVWPAVVFTLVTLLDDPVAADLAVDAVEPVATAWDTAARIGLADRRLYLAANRCLAIAARRVPTELIGAMQRLVDHVDRGVCPADDFSDRVIAGGIASAVTGMMHGAS",
   "taxon_ids": [
    "83332"
   ],
   "accessions": [
    "A0A1R3Y503",
    "A0ABX2VMI8",
    "A0ACD3TZ38",
    "A5U929",
    "P9WPK6",
    "P9WPK7"
   ],
   "assertions": [
    {
     "assertion_id": "97f334e055048d81c04ba0a6a113c16782448fe9ff0db2dc54f20cfd528f2794:GO:0005524",
     "go_id": "GO:0005524",
     "go_name": "ATP binding",
     "go_definition": "\"Binding to ATP, adenosine 5'-triphosphate, a universally important coenzyme and enzyme regulator.\" [ISBN:0198506732]",
     "supporting_assertions": [
      {
       "accession": "A5U929",
       "symbol": "gshA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:419947"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0ABX2VMI8",
       "symbol": "egtA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:1844474"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "P9WPK6",
       "symbol": "egtA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:83331"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0ACD3TZ38",
       "symbol": "egtA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:1305738"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "P9WPK7",
       "symbol": "egtA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:83332"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      },
      {
       "accession": "A0A1R3Y503",
       "symbol": "gshA",
       "relation": "enables",
       "go_id": "GO:0005524",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000184798"
       ],
       "object_name": "Glutamate--cysteine ligase EgtA",
       "taxon": [
        "taxon:233413"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "ae7a8231780df578d57ef36d2619bf4ff6c8d44318204a57f3a04de697f06eed",
   "sequence": "MPELPEVEAARRAIEENCLGKKIKRVIIADDNKVIHGISPSDFQTSILGKTIISARRKGKNLWLELDSPPFPSFQFGMAGAIYIKGVAVTKYKRSAVKDSEEWPSKYSKFFVELDDGLELSFTDKRRFAKVRLLANPTSVSPISELGPDALLEPMTVDEFAESLAKKKITIKPLLLDQGYISGIGNWIADEVLYQARIHPLQTASSLSKEQCEALHTSIKEVIEKAVEVDADSSQFPSYWIFHNREKKPGKAFVDGKKIDFITAGGRTTAYVPELQKLYGKDAEKAAKVRPAKRGVKPKEDDGDGEEDEQETEKEDESAKSKKGQKPRGGRGKKPASKTKTEESDDDGDDSEAEEEVVKPKGRGTKPAIKRKSEEKATSQAGKKPKGRKS",
   "taxon_ids": [
    "3702"
   ],
   "accessions": [
    "O80358"
   ],
   "assertions": [
    {
     "assertion_id": "ae7a8231780df578d57ef36d2619bf4ff6c8d44318204a57f3a04de697f06eed:GO:0003684",
     "go_id": "GO:0003684",
     "go_name": "damaged DNA binding",
     "go_definition": "\"Binding to damaged DNA.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "O80358",
       "symbol": "FPG1",
       "relation": "enables",
       "go_id": "GO:0003684",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR015886"
       ],
       "object_name": "Formamidopyrimidine-DNA glycosylase",
       "taxon": [
        "taxon:3702"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR015886/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "1c71b14ddab17afc222e0d44923249d7fceda757299283f7aedb508a707bac82",
   "sequence": "MIMEIPAIKALSRYAQWVIWKKERDTKIPYNPNNGKKASSTDPLAWGDIDEAQAGLVRYGANGLGFVLTKSDPFVFIDLDHVLDENKRVKCEWARQLLKEIKSYTEISPSGDGLHVVVSGKLPDYIKHKTKFDDGSALEVYESGRYMTITGEVFDGRDDIKELDLSILGEFAEHKIETKNAPVQIESATTLDDEAIIDLMKRKGQWPDAPKDGDDWSSLDMSFANRLAFWCGKDIERMDRIFRQSPLMRQKWDRPTAGSTYGRITLKKACDFVDSVYDPALRNESDCPFEPYNEEGGPRNDKEEKDPLWLYKVLLTKGIEVWFDIKLEKYGIKRNNRVDYIAKSSLQQIVFEIIGKTPKNIAVPTYIGAYEPSKPEKWEEEGIKYINLFKPTPLMKVKPVKEMPEIVKNLLLNLFDYDAKSMGLFINWLAFIYQYKERTGVAWIFMGKQGTGKGLLVDLLKKIFEEHMSSNITDANLDSQFNPYLYNKLIVHLNEVSADNRKSRMLVKNRLKTWITDETLYINRKNMKEVEIKNFCNFIINSNETIPVDIEDSDRRFNVIECNNVLKEQEWWTTESYQEILNNAEGFAKYLAGIKVDRSKVNEVVMSEKKKAIVETTESVLKQIAKALTDRDIEWFLDNGLEGVVEKNIVNDFQWEELQEAITTGVIPNKYLMIIVEQILGDSKTITWIKRNIITPYQVGETTVVKMAGKPIRAIVVG",
   "taxon_ids": [
    "1230469"
   ],
   "accessions": [
    "M5AAG8"
   ],
   "assertions": [
    {
     "assertion_id": "1c71b14ddab17afc222e0d44923249d7fceda757299283f7aedb508a707bac82:GO:0003678",
     "go_id": "GO:0003678",
     "go_name": "DNA helicase activity",
     "go_definition": "\"Unwinding of a DNA helix, driven by ATP hydrolysis.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "M5AAG8",
       "symbol": "28",
       "relation": "enables",
       "go_id": "GO:0003678",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:3.6.4.12"
       ],
       "object_name": "DNA Primase-polymerase",
       "taxon": [
        "taxon:1230469"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "a5accf714119c3a78523d6e5b48f1816d902b31f8107af60f41fced1ba43e758",
   "sequence": "MSSSIPRVYSLGNSAMTYLLALRIAQLPSQPKVPSVVLLLNDQKKLNRFLNNDSKIIVKSSNNNKETYHRQFMASCVPPILSNGEIAPIENLIVSDPSSKFITAQLSKYNKSLRPETNILFLNPSLNLLEHLHRYRWRFDEARPNLFMGFTTPVDVGTIHQEFQLSLKVKGRIQFHIAKIDGFPRMSSTGKSASLSLRGDRQKNEKENNAFYKLFREISRLRSGIGSDLVSFDLHVHGFQDLFFTELEKLILESCTEPLLAVYDCVYKKELLKIPGAQDIIKKLISEQLSIIDRSYPSLNTNPNYSVIFDKERIFSLVMRDLEVNGHKRAKLAQSLNQLNQTNINELNGFFVSLGKYKKCNCKWNDILLTLIKGKQFITKQKALDYHYL",
   "taxon_ids": [
    "559292"
   ],
   "accessions": [
    "B3LG86",
    "H0GDL5",
    "P14905"
   ],
   "assertions": [
    {
     "assertion_id": "a5accf714119c3a78523d6e5b48f1816d902b31f8107af60f41fced1ba43e758:GO:0050661",
     "go_id": "GO:0050661",
     "go_name": "NADP binding",
     "go_definition": "\"Binding to nicotinamide-adenine dinucleotide phosphate, a coenzyme involved in many redox and biosynthetic reactions; binding may be to either the oxidized form, NADP+, or the reduced form, NADPH.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "H0GDL5",
       "symbol": "VIN7_0912",
       "relation": "enables",
       "go_id": "GO:0050661",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN000496353"
       ],
       "object_name": "Ketopantoate reductase C-terminal domain-containing protein",
       "taxon": [
        "taxon:1095631"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      },
      {
       "accession": "B3LG86",
       "symbol": "SCRG_00323",
       "relation": "enables",
       "go_id": "GO:0050661",
       "references": [
        "GO_REF:0000118"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "PANTHER:PTN000496353"
       ],
       "object_name": "Ketopantoate reductase C-terminal domain-containing protein",
       "taxon": [
        "taxon:285006"
       ],
       "annotation_date": "20260727",
       "assigned_by": "TreeGrafter",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000118.html"
       ]
      },
      {
       "accession": "P14905",
       "symbol": "CBS2",
       "relation": "enables",
       "go_id": "GO:0050661",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "PANTHER:PTN002464833",
        "UniProtKB:P0A9J4"
       ],
       "object_name": "Cytochrome B translational activator CBS2",
       "taxon": [
        "taxon:559292"
       ],
       "annotation_date": "20180502",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/P0A9J4/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "e2ee6aac82b25888a5428eaa5daee45e72b9aa61de6469bcaa00f15695672b4c",
   "sequence": "MATARGGSGPDPGSRGLLLLSFSVVLAGLCGGNSVERKIYIPLNKTAPCVRLLNATHQIGCQSSISGDTGVIHVVEKEDDLKWVLTDGPNPPYMVLLEGKLFTRDIMEKLKGETSRIAGLAVTLAKPNSTSSFSPSVQCPNDGFGIYSNSYGPEFAHCKKTLWNELGNGLAYDDFSFPIFLLEDENETKVIKQCYQDHNLGQNGSAPSFPLCAMQLFSHMHAVISTATCMRRSFIQSTFSINPEIVCDPLSDYNVWSMLKPINTSGGLEPDVRVVVAATRLDSRSFFWNVAPGAESAVASFVTQLAAAEALHKAPDVTTLPRNVMFVFFQGETFDYIGSSRMVYDMENGKFPVRLENIDSFVELGQVALRTSLELWMHTDPMSQKNESVKNQVEDLLVTLEQSGADTPQVVLSRLVQSQALPPSSLQRFLRARNISGVVLADHSGSFHNRYYQSIYDTAENINVTYPESQSPEEDLNFVTDTAKALADVATVLARALYKLAGGTNFNNSIQADPQTVTRLLYGFLVRANNSWFQSILRHDLRSYLDDGPLQHYIAVSSPTNTTYVVQYALANLTGKVTNLTQEQCQDPSKVPNESKDLYEYSWVQGPWNSNKTERLPRCVRSTVRLARALSPAFELSQWSSTEYSTWAESRWKDIQARIFLIASKELEFITLIVGFSILVFSLIVTYCINAKADVLFVAPREPGAVSY",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "Q8CGU6"
   ],
   "assertions": [
    {
     "assertion_id": "e2ee6aac82b25888a5428eaa5daee45e72b9aa61de6469bcaa00f15695672b4c:GO:0030674",
     "go_id": "GO:0030674",
     "go_name": "protein-macromolecule adaptor activity",
     "go_definition": "\"An adaptor activity that brings together two or more macromolecules in contact, permitting those molecules to function in a coordinated way. The adaptor can bring together two proteins, or a protein and another macromolecule such as a lipid or a nucleic acid.\" [GOC:bf, GOC:mah, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "Q8CGU6",
       "symbol": "Ncstn",
       "relation": "enables",
       "go_id": "GO:0030674",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q92542"
       ],
       "object_name": "Nicastrin",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20201119",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/Q92542/entry"
       ]
      },
      {
       "accession": "Q8CGU6",
       "symbol": "Ncstn",
       "relation": "enables",
       "go_id": "GO:0030674",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q92542",
        "ensembl:ENSP00000294785"
       ],
       "object_name": "Nicastrin",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q92542/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "56a5c62fee4fbf845b49154edd33e50306e89013b89a6a31850315febe36204d",
   "sequence": "MFTHLTRAFRKMNNLVNRSFIDVHRVVAELSYPEFEEDVKNPESSIYRTPISLFQNKDIVTIVGDYILSPKTDSFQVLYPIKKVIEHFPVIFHCTHNNAPLWVHLLDERHHRLLQSLLTYEIVNAKYRGIVVIPYYRRPINYQTGKSLLMSKLASVKVLDILMRCGSYKFISLMCMINKKNNTNFLHCCASKWGEVGSKMMLHIAEMFFANPTTSQHLSDASSFPDAAAEDDKGKTPAHLAIQEDNADALLFLISLYGAPWFQDNNSYMKSALELKSNKCVKVLSFAADKYEILPNINNNQLEPDTMCGVCATSVEEDENEGKTTSLSWYQMNCKHYIHCECLMGMCAAAGNVQCPMCREDVGDEVLERCPPTIFRWLKLAERSEHNRVLFEAKKQEFYKQMEAMKPPRVVVPPRRTFLTPARRGERAIRIAREIATNAIAEATAQGDVNSYFPVLIDGSGEEYEEEGEEFFNSEEEALAFGRPFLEDEEEARQIQMRQFAELSRRGVSVNIINNDNPHRHISTVNIVQPVYGVEKSPAASFIYNMLKNDVFESIRSRDTRVGGERVPVMNLSNDKRALFHAASSMLCDFATETNSQIVGLDFQAVYDPHHISNYIETFGSPLHAYPGAVTFLDGAQDYYAESIRYDNDIVSFSEMASELHITEALDVFEGSLLSPLFKKIRTGKSYSNWNDHLRRRNYARDIAEEFVRVCENSLASREHPPVHVHPFRDGAIPILIEYIVDFIHHCITWSMQVNALHCMRKYIEHENTNVHLLNLRPTDERVEVLRVSQLRWSRLFNEQYNTRMSLSTKRLSLMKIFNHDLGVSKFGVYKLLDIIEMYCFTLI",
   "taxon_ids": [
    "654913"
   ],
   "accessions": [
    "A0A0S2E6R1",
    "A0ACM8Q1H8",
    "Q77F08",
    "Q77J49",
    "Q91LD2"
   ],
   "assertions": [
    {
     "assertion_id": "56a5c62fee4fbf845b49154edd33e50306e89013b89a6a31850315febe36204d:GO:0000166",
     "go_id": "GO:0000166",
     "go_name": "nucleotide binding",
     "go_definition": "\"Binding to a nucleotide, any compound consisting of a nucleoside that is esterified with (ortho)phosphate or an oligophosphate at any hydroxyl group on the ribose or deoxyribose.\" [GOC:mah, ISBN:0198547684]",
     "supporting_assertions": [
      {
       "accession": "Q77J49",
       "symbol": "WSV222",
       "relation": "enables",
       "go_id": "GO:0000166",
       "references": [
        "GO_REF:0000043"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB-KW:KW-0547"
       ],
       "object_name": "RING finger containing E3 ubiquitin-protein ligase WSV222",
       "taxon": [
        "taxon:654913"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000043.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "e2f1ad690491b11f8639df861eec336bdcc78d0777504d6ec59e70bdc379f237",
   "sequence": "MLLSVPLLLGLLGLAAADPAIYFKEQFLDGDAWTNRWVESKHKSDFGKFVLSSGKFYGDQEKDKGLQTSQDARFYALSARFEPFSNKGQTLVVQFTVKHEQNIDCGGGYVKLFPGGLDQKDMHGDSEYNIMFGPDICGPGTKKVHVIFNYKGKNVLINKDIRCKDDEFTHLYTLIVRPDNTYEVKIDNSQVESGSLEDDWDFLPPKKIKDPDAAKPEDWDERAKIDDPTDSKPEDWDKPEHIPDPDAKKPEDWDEEMDGEWEPPVIQNPEYKGEWKPRQIDNPDYKGTWIHPEIDNPEYSPDANIYAYDSFAVLGLDLWQVKSGTIFDNFLITNDEAYAEEFGNETWGVTKAAEKQMKDKQDEEQRLKEEEEDKKRKEEEEAEDKEDEDDRDEDEDEEDEKEEDEEDATGQAKDEL",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "P18418"
   ],
   "assertions": [
    {
     "assertion_id": "e2f1ad690491b11f8639df861eec336bdcc78d0777504d6ec59e70bdc379f237:GO:0030674",
     "go_id": "GO:0030674",
     "go_name": "protein-macromolecule adaptor activity",
     "go_definition": "\"An adaptor activity that brings together two or more macromolecules in contact, permitting those molecules to function in a coordinated way. The adaptor can bring together two proteins, or a protein and another macromolecule such as a lipid or a nucleic acid.\" [GOC:bf, GOC:mah, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "P18418",
       "symbol": "Calr",
       "relation": "enables",
       "go_id": "GO:0030674",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:P27797"
       ],
       "object_name": "Calreticulin",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260312",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/P27797/entry"
       ]
      },
      {
       "accession": "P18418",
       "symbol": "Calr",
       "relation": "enables",
       "go_id": "GO:0030674",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:P27797",
        "ensembl:ENSP00000320866"
       ],
       "object_name": "Calreticulin",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/P27797/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "689eec4d6ec2210cc2b7956699a74e51429ef6413da127acf99a5d721e551e8a",
   "sequence": "MRRSSVTNRQSLLPLRVQDANRMGLTTPQNKDRQGFGKLSMSKPHSGTSEKKTSFFGKRASNGTARTSQYGAFGGTEKIKDPRPLHDKAFIQQCIRQLCEFLNENGYSQALTVKSLQGPSTKDFLKIFAFIYTFICPNYENPESKFEEEIPRIFKELGYPFALSKSSMYTVGAPHTWPQIVAALVWLIDCVKLCCVLRSENPVFEDPPMGEQSENGIDFNQLFLDYTVRCYDQFMEGRDTFEEYDTDVCIRLKEAYHVDESNLEALQQESRRLMEEIERLEKEKENEPDRLASMRKLKASLQADIQKYQNYLTEIESHSTLLDQRVSSISEELEATELESRAIQQENLRLKNILDNQKYSVADIERIKYEETELQQTIAKLTKDLDEDKQQLWHEELKYAKMKESVETQLSEFHKIGRKVRLIPPTAEFANGYDFQIQCNLDSEQSSLMHYRNKINVPLVEILSQSEGHIASSTNKKMGVEDMVEQLNTLIGEKKNEVKVQKDEAQKLEEIYQQKVEESEEHEKRWISEIESLEKHRQLLESGVNKSLDEAMKDLQKEQQELQLVEHRTEEEMRQVENKLVRVVHAVANHVAVIEKHLEEKRLKVEREYEEFMKEDLLLDLRELLEKYKEKARVLDSL",
   "taxon_ids": [
    "8355"
   ],
   "accessions": [
    "Q8AWF5"
   ],
   "assertions": [
    {
     "assertion_id": "689eec4d6ec2210cc2b7956699a74e51429ef6413da127acf99a5d721e551e8a:GO:0008017",
     "go_id": "GO:0008017",
     "go_name": "microtubule binding",
     "go_definition": "\"Binding to a microtubule, a filament composed of tubulin monomers.\" [GOC:krc]",
     "supporting_assertions": [
      {
       "accession": "Q8AWF5",
       "symbol": "ndc80",
       "relation": "enables",
       "go_id": "GO:0008017",
       "references": [
        "GO_REF:0000024"
       ],
       "evidence_code": "ISS",
       "with_from": [
        "UniProtKB:O14777"
       ],
       "object_name": "Kinetochore protein NDC80 homolog",
       "taxon": [
        "taxon:8355"
       ],
       "annotation_date": "20230809",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000024.html",
        "https://www.uniprot.org/uniprotkb/O14777/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "0c7464842c968e3865a429a57d0722e80c2736068de105c0f014a046c6cd98c6",
   "sequence": "MRIFYYLHFLCYVTFILPATCTLVNADRCTKRYGRCKRDCLESEKQIDICSLPRKICCTEKLYEEDDMF",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "Q30KQ6"
   ],
   "assertions": [
    {
     "assertion_id": "0c7464842c968e3865a429a57d0722e80c2736068de105c0f014a046c6cd98c6:GO:0031731",
     "go_id": "GO:0031731",
     "go_name": "CCR6 chemokine receptor binding",
     "go_definition": "\"Binding to a CCR6 chemokine receptor.\" [GOC:mah, GOC:nln]",
     "supporting_assertions": [
      {
       "accession": "Q30KQ6",
       "symbol": "DEFB114",
       "relation": "enables",
       "go_id": "GO:0031731",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "MGI:MGI:1927667",
        "PANTHER:PTN000464032",
        "UniProtKB:O15263"
       ],
       "object_name": "Beta-defensin 114",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20170427",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html",
        "https://www.uniprot.org/uniprotkb/O15263/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "4f8761cb76fd9e542d25a37c438a64f2c8c6d8d75884a989504f9d2c3d9541ed",
   "sequence": "MGDKINFIFGIHNHQPLGNFGWVFEEAYEKCYWPFLETLEEYPNMKVAIHTSGPLIEWLQDNRPEYIDLLRSLVKRGQVEIVVAGFYEPVLASIPKEDRIEQIRLMKEWAKSIGFDARGVWLTERVWQPELVKTLKESGIDYVIVDDYHFMSAGLSKEELYWPYYTEDGGEVIAVFPIDEKLRYLIPFRPVDKVLEYLHSLIDGDESKVAVFHDDGEKFGIWPGTYEWVYEKGWLREFFDRISSDEKINLMLYTEYLEKYKPRGLVYLPIASYFEMSEWSLPAKQARLFVEFVNELKVKGIFEKYRVFVRGGIWKNFFYKYPESNYMHKRMLMVSKLVRNNPEARKYLLRAQCNDAYWHGLFGGVYLPHLRRAIWNNLIKANSYVSLGKVIRDIDYDGFEEVLIENDNFYAVFKPSYGGSLVEFSSKNRLVNYVDVLARRWEHYHGYVESQFDGVASIHELEKKIPDEIRKEVAYDKYRRFMLQDHVVPLGTTLEDFMFSRQQEIGEFPRVPYSYELLDGGIRLKREHLGIEVEKTVKLVNDGFEVEYIVNNKTGNPVLFAVELNVAVQSIMESPGVLRGKEIVVDDKYAVGKFALKFEDEMEVWKYPVKTLSQSESGWDLIQQGVSYIVPIRLEDKIRFKLKFEEASG",
   "taxon_ids": [],
   "accessions": [
    "P49067"
   ],
   "assertions": [
    {
     "assertion_id": "4f8761cb76fd9e542d25a37c438a64f2c8c6d8d75884a989504f9d2c3d9541ed:GO:0004556",
     "go_id": "GO:0004556",
     "go_name": "alpha-amylase activity",
     "go_definition": "\"Catalysis of the endohydrolysis of (1->4)-alpha-D-glucosidic linkages in polysaccharides containing three or more alpha-(1->4)-linked D-glucose units.\" [PMID:12527308]",
     "supporting_assertions": [
      {
       "accession": "P49067",
       "symbol": "amyA",
       "relation": "enables",
       "go_id": "GO:0004556",
       "references": [
        "GO_REF:0000003"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "EC:3.2.1.1"
       ],
       "object_name": "Alpha-amylase",
       "taxon": [
        "taxon:186497"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000003.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "3038428ff8d5bfc09438641103de1e40984d86cac380264b1748d5c9c345b61a",
   "sequence": "MEVRKINTHNREKKLDDIVYDDFISTLRQIKEGNHQLREEFISEYKPFILKVTSNATGKYIDTRNSDEFSIALSAFNEAIDKFDIEKGYNFFLFSEQVIRRRLIDYSRSNKDDKEYPFSFFDDEYFYNNEKLLSKSYIGFEDIEAREDIEELKKKLQEFGITFLDLVLNVPKHRDSRQLCIRLAKMLAEDEQMYNALMKNKNIPRNELKKKAKVHGRTIGNNRKYIIALCLIFRSNLNLSKRYLEYYTMDGESDLI",
   "taxon_ids": [
    "203119"
   ],
   "accessions": [
    "A3DBH0"
   ],
   "assertions": [
    {
     "assertion_id": "3038428ff8d5bfc09438641103de1e40984d86cac380264b1748d5c9c345b61a:GO:0003677",
     "go_id": "GO:0003677",
     "go_name": "DNA binding",
     "go_definition": "\"Any molecular function by which a gene product interacts selectively and non-covalently with DNA (deoxyribonucleic acid).\" [GOC:dph, GOC:jl, GOC:tb, GOC:vw]",
     "supporting_assertions": [
      {
       "accession": "A3DBH0",
       "symbol": "sigI1",
       "relation": "enables",
       "go_id": "GO:0003677",
       "references": [
        "GO_REF:0000104"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniRule:UR000336380"
       ],
       "object_name": "RNA polymerase sigma factor SigI1",
       "taxon": [
        "taxon:203119"
       ],
       "annotation_date": "20260727",
       "assigned_by": "UniProt",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000104.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "a20347dfb24db3b0e3878f27902d11ec84e5754746c4c7ce676ceb38de7fe1c2",
   "sequence": "MASAPLPGPPASAGGDGPNLNNNNNNNNNHSVRKCGYLRKQKHGHKRFFVLRGPGTGGEEAAAAGGSPPQPPRLEYYESEKKWKSKAGAPKRVIALDCCLNINKRADAKHKYLIALYTKDEYFAVAAENEQEQEGWYRALTDLVSEGRSGDGGSGTTGGSCSASLPGALGGSAGAAGCDDNYGLVTPATAVYREVWQVNLKPKGLGQSKNLTGVYRLCLSARTIGFVKLNCEQPSVTLQLMNIRRCGHSDSFFFIEVGRSAVTGPGELWMQADDSVVAQNIHETILEAMKALKELFEFRPRSKSQSSGSSATHPISVPGARRHHHLVNLPPSQTGLVRRSRTDSLAATPPAAKCTSCRVRTASEGDGGAAGGAGTAGGRPMSVAGSPLSPGPVRAPLSRSHTLSAGCGGRPSKVALAPAGGALQHSRSMSMPVAHSPPAATSPGSLSSSSGHGSGSYPLPPGSHPHLPHPLHHPQCQRPSSGSASASGSPSDPGFMSLDEYGSSPGDLRAFSSHRSNTPESIAETPPARDGSGGELYGYMSMDRPLSHCGRPYRRVSGDGAQDLDRGLRKRTYSLTTPARQRQVPQPSSASLDEYTLMRATFSGSSGRLCPSLPASSPKVAYNPYPEDYGDIEIGSHKSSSSNLGADDGYMPMTPGAALRSGGPNSCKSDDYMPMSPTSVSAPKQILQPRSAAALPPSGAAVPAPPSGAGRTFPVNGGGYKASSPAESSPEDSGYMRMWCGSKLSMENPDPKLLPNGDYLNMSPSEAGTAGTPPDFFSAALRPGGEALKGVPGHCYSSLPRSYKAPCTCGGGDNDQYVLMSSPVGRILEEERLEPQATPGAGTFGAAGGSHTQPHHSAVPSSMRPSGIVGRPEGFLGQRCRAVRPTRLSLEGLQTLPSMQEYPLPTEPKSPGEYINIDFGEGGTRLSPPAPPLLASAASSSSLLSASSPASSLGSGTPGTSSDSRQRSPLSDYMNLDFSSPKSPKPSTRSGDTVGSIDGLLSPEASSPYPPLPPRPSASPSSLQQPLPPAPGDLYRLPPATAATSQGPTAGSSMSSEPGDNGDYTEMAFGVAATPPQPIAAPPKPEGARVTSPTSGLKRLSLMDQVSGVEAFLQVSQPPDPHRGAKVIRADPQGGRRRHSSETFSSTTTVTPVSPSFAHNSKRHNSASVENVSLRKSSEGNSILGGSDEPSTSPGQAQPSAGVPPAPQARPWNPGQPGALIGCPGGSSSPMRRETSVGFQNGLNYIAIDVRGEQGSLAQSQPQHPQPGDKNSWGRTRSLGGLLGTVGGSGTSGVCGGPGTGALPSASTYASIDFLSHHLKEATVVKE",
   "taxon_ids": [
    "10116"
   ],
   "accessions": [
    "F1MAL5"
   ],
   "assertions": [
    {
     "assertion_id": "a20347dfb24db3b0e3878f27902d11ec84e5754746c4c7ce676ceb38de7fe1c2:GO:0035591",
     "go_id": "GO:0035591",
     "go_name": "signaling adaptor activity",
     "go_definition": "\"The binding activity of a molecule that brings together two or more molecules in a signaling pathway, permitting those molecules to function in a coordinated way. Adaptor molecules themselves do not have catalytic activity.\" [GOC:bf, PMID:19104498]",
     "supporting_assertions": [
      {
       "accession": "F1MAL5",
       "symbol": "Irs2",
       "relation": "enables",
       "go_id": "GO:0035591",
       "references": [
        "GO_REF:0000121"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q9Y4H2"
       ],
       "object_name": "Insulin receptor substrate 2",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20240731",
       "assigned_by": "RGD",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000121.html",
        "https://www.uniprot.org/uniprotkb/Q9Y4H2/entry"
       ]
      },
      {
       "accession": "F1MAL5",
       "symbol": "Irs2",
       "relation": "enables",
       "go_id": "GO:0035591",
       "references": [
        "GO_REF:0000107"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "UniProtKB:Q9Y4H2",
        "ensembl:ENSP00000365016"
       ],
       "object_name": "Insulin receptor substrate 2",
       "taxon": [
        "taxon:10116"
       ],
       "annotation_date": "20260724",
       "assigned_by": "Ensembl",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000107.html",
        "https://www.uniprot.org/uniprotkb/Q9Y4H2/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "367701a601ecb0202fd500db8dbac445812f1cc03e500ca98774209b70b0ce69",
   "sequence": "MGPSRDACDRCHTMKTRCQRTPQTRECVRCNRLGISCTYSPPGRTGRPLGARKGGSEKGKKEAAARSRKGGRISPMVKTHMSTSLGMPTPGLDHSLQMMPDGNSGVNPSLLGDFDLFSASAYQAELDYLSFPFMDSFLAESAGDLTFAFASAMDSQPSYSSGSSSPPGPSSPRTPATPHQQPFQPSSSHDERAWSVNNTDSFSDPGSPSSLFTENIFSPSFIPKPPRNNQYFPSDERSDNGEGDGEVLVRLLDLQARISKLVKHLHGGRGSSIDCEEVLNAGKSLISLLDSVGAPSRFFSSSSSSPFSDSRMSEPASPMFSSSSRSTSSSTSSALLSSENTIPSNVLILSLSSCYATLLHAYDLLADRLHHHHQQQHSHSHCQKYVSGIQTPLLSPASSISSFSPIEQQQGYFFGNGLKHQAQQPLDGGNNAQLHAMAQMVQKLKNTLQCFMGRLGQGSSGANGALGAKKKKKQNPNVHLDMESCCGNASNDDDLLMSLSKPALLEMGWGEEELWQGMNW",
   "taxon_ids": [],
   "accessions": [
    "Q1K836"
   ],
   "assertions": [
    {
     "assertion_id": "367701a601ecb0202fd500db8dbac445812f1cc03e500ca98774209b70b0ce69:GO:0000981",
     "go_id": "GO:0000981",
     "go_name": "DNA-binding transcription factor activity, RNA polymerase II-specific",
     "go_definition": "\"A DNA-binding transcription factor activity that modulates the transcription of specific gene sets transcribed by RNA polymerase II.\" [GOC:txnOH-2018]",
     "supporting_assertions": [
      {
       "accession": "Q1K836",
       "symbol": "NCU01209",
       "relation": "enables",
       "go_id": "GO:0000981",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR001138",
        "InterPro:IPR036864"
       ],
       "object_name": "Zn(2)-C6 fungal-type domain-containing protein",
       "taxon": [
        "taxon:367110"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001138/",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR036864/"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "73135a5306ecb12c0e7f2802a8a39520dc6e427d65bfd4b12e29831bb1faed97",
   "sequence": "MDITVSELMSNFMDSPLVVWVKTFGPLGFSSEDKLSMFMDLVDGVFLHKIMTHIDPSPMNQRVNKQVNNDVNLRIQNLNTVIRHIKNYYQEHLQQLIVMNLPNVLAIAKDPLSGKSMEEMKRMLLLILGCAVQCDRKEEIIEKIKLLDIETQAAIVTHIQEVTHNQENVLDLQWMEVAEIPAEQLDPLSRTMAFHLRKLIDERDESAELVIELTQERDYLQSQQPSGLLGFPSPERTSLSPITLLSKEDRQHLAVELADTKAKLRRSRQELEEKTEQLIDAKNEIERLDSDIQKLKQENTQLLAEARSVRAYRDEVDSLRERAGKVDRLETELSRFKEKLNDVHFYKTRIEELREDNLTLLETKSMLEEQLTGARGRCDKLHELEKENLQLRSKLHDIEIDRDSDKKRLEELLEENMLLEISQKQSMNESAHLGWELEQLAKNNEVNEARKSFVFELNESASSRLLKLEKENQCLQSTIQELREASINMEEGQLHSLELEKENQSLSKKLERLQSQLDQEKQTTQDMENLGEELIKEKQRMEKTLETIQAEKDRQISELEQEKEHLTQAVSSLRKRAQANSEARVREVETENRILHQTISETGGKLARLEAEKRQVTKELESLRERGERCEELEREVPRLERVREQLQREAAALKIGSERAEALERENATLEQDNRRLKKLADTAQNATLRLAVLEKDHQQLEEENLEQRRALETLRPAAARLAQLQQEHAELEREHEEMCRTMEELRSQAKRSERLEKSCGSLSLENQRLQQTLENSSTKMQGLESELRQNEAEMKDLQRELEGLRQKVTWAETLEKENRGMEQELSQLEKEKKQLEKEARRFRQQLEVKEAALEENCLRLASMEKEGTALSKELGRVKEAAGRLKELERENKDLQKQATMDKKTLATLREELVNEKLRVQQQCNELEKLSHELEKIGLNREKLLQEEHSCEDNKYKILETKIESALKKTLELREEKIQSLESRLEESSSLNQQLRTELTTVKKNLEALKQRHEEEAAHSEISQQTLGQTRSLPDKEKWEMEQREATAELLKLKDRLIDVEKNVRQRHVSIDIHRVIFSIVICFCDSLQNAALQTEKYLLKDQLKQIDSQNAQLNAQTLALQKQAASLQEHNTSLHKETAKLQVENSTLSSQSSSLMAQYGALQAQLQTLESEAESLQKQREEASAARDRVTQDHERLLGVHERQASEYEQLIAQHAALKASQRALEQENRTLENKYMVLLKQKDAMEALEESLQRDRESLGEEIRKNTLILGENRSLREEVDRVSHMHTQLRQEYDSLQLQTKELKTSLNESQLELNRWQARYDQLKEQHQGLDISMTKLDNHCELLTRLKGNLEEENHHLLSQIQMLSQQNQTLLERTMESKELYHEEQKQYIDKLNSLRRQKEKLEEKIMDQYKFYDPTPKKSRQWVGAKAIAKFIKPKKESSRERPDAPRERIRSAPDIPLPEIPTCIDCPESAPPPPPPPLPPRQSRPSLDSMNSQSVEENHVQSPTLSSPALNGRVLNESGGSRSRDGYRSIGGGSESMNGYEELLRWRSREPGGATCSTPLSRNSHNAPGFTSSSSLRPGRRPKGLVSEEDLRHHSPDAGFGSGVHGNTGHRPSSAEFSRNTSSSNSPVSSKGSLDCLQGRSASLSSDDVVGLAHEGSRLSQSSLLPRSSTLPCDSPSASRPSQRPASRRPSSPGSEMVTLEEFLQESNALSPPTVQTGSREDLMTDYFTRSTRPVPLRDGAKTPTNYVTPTVKTTPPELDARTPKPGHSVKPSVRLTDTSTPPSHSQTLPNRGAGLRPSALQQSSPRGSVGGSASLSRTFSLASADLLRSNGPDSYRTEAASPNQNDVVMRRPGAVARERPMSARVTGSSPLPGDPGHISVDPRRLSLAQPRDEFSLVSPPPLHSSSMSLQAEREYVGSGSSRAGAARSGSAQPRGAPHRGEVAMVTPVRAVPALRLNDLEEEPQEQREAESPLLKKADTTNLSYASKEQPTSKPASPDPNNDPQTVWYEYGCV",
   "taxon_ids": [
    "7955"
   ],
   "accessions": [
    "A0A2R8QCI3"
   ],
   "assertions": [
    {
     "assertion_id": "73135a5306ecb12c0e7f2802a8a39520dc6e427d65bfd4b12e29831bb1faed97:GO:0051959",
     "go_id": "GO:0051959",
     "go_name": "dynein light intermediate chain binding",
     "go_definition": "\"Binding to a light intermediate chain of the dynein complex.\" [GOC:bf]",
     "supporting_assertions": [
      {
       "accession": "A0A2R8QCI3",
       "symbol": "ccdc88c",
       "relation": "enables",
       "go_id": "GO:0051959",
       "references": [
        "GO_REF:0000033"
       ],
       "evidence_code": "IBA",
       "with_from": [
        "PANTHER:PTN000436622",
        "WB:WBGene00006997"
       ],
       "object_name": "Protein Daple",
       "taxon": [
        "taxon:7955"
       ],
       "annotation_date": "20170228",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000033.html"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "6c1efe2ee9c429e4832d89e49c0067c32db4c7400beed36961e6ba1035e4d66c",
   "sequence": "MMNEDAAQKSDSGEKFNGSSQRRKRPKKSDSNASFLRAARAGNLDKVVEYLKGGIDINTCNQNGLNALHLAAKEGHVGLVQELLGRGSSVDSATKKGNTALHIASLAGQAEVVKVLVKEGANINAQSQNGFTPLYMAAQENHIDVVKYLLENGANQSTATEDGFTPLAVALQQGHNQAVAILLENDTKGKVRLPALHIAARKDDTKSAALLLQNDHNADVQSKMMVNRTTESGFTPLHIAAHYGNVNVATLLLNRGAAVDFTARNGITPLHVASKRGNTNMVKLLLDRGGQIDAKTRDGLTPLHCAARSGHDQVVELLLERKAPLLARTKNGLSPLHMAAQGDHVECVKHLLQYKAPVDDVTLDYLTALHVAAHCGHYRVTKLLLDKRANPNARALNGFTPLHIACKKNRIKVMELLVKYGASIQAITESGLTPIHVAAFMGHLNIVLLLLQNGASPDVTNIRGETALHMAARAGQVEVVRCLLRNGALVDARAREEQTPLHIASRLGKTEIVQLLLQHMAHPDAATTNGYTPLHISAREGQVDVASVLLEAGAAHSLATKKGFTPLHVAAKYGSLDVAKLLLQRRAAADSAGKNGLTPLHVAAHYDNQKVALLLLEKGASPHATAKNGYTPLHIAAKKNQMQIASTLLNYGAETNTVTKQGVTPLHLASQEGHTDMVTLLLDKGANIHMSTKSGLTSLHLAAQEDKVNVADILTKHGADRDAYTKLGYTPLIVACHYGNVKMVNFLLKQGANVNAKTKNGYTPLHQAAQQGHTHIINVLLQHGAKPNATTANGNTALAIAKRLGYISVVDTLKVVTEEVTTTTTTITEKHKLNVPETMTEVLDVSDEEGDDTVTGDGGEYLRPEDLKELGDDSLPSSQFLDGMNYLRYSLEGGRSDSLRSFSSDRSHTLSHASYLRDSAMIDDTVVIPSHQVSALAKEAERNSYRLSWGTENLDNVALSSSPIHSGFLVSFMVDARGGAMRGCRHNGLRIIIPPRKCTAPTRVTCRLVKRHRLATMPPMVEGEGLASRLIEVGPSGAQFLGPVIVEIPHFAALRGKERELVVLRSENGDSWKEHFCDYTEDELNEILNGMDEVLDSPEDLEKKRICRIITRDFPQYFAVVSRIKQDSNLIGPEGGVLSSTVVSQVQAVFPEGALTKRIRVGLQAQPMHSELVKKILGNKATFSPIVTLEPRRRKFHKPITMTIPVPKASSDVMLNGFGGDAPTLRLLCSITGGTTPAQWEDITGTTPLTFVNECVSFTTNVSARFWLIDCRQIQESVAFASQVYREIICVPYMAKFVVFAKSHDPIEARLRCFCMTDDKVDKTLEQQENFSEVARSRDVEVLEGKPIYVDCFGNLVPLTKSGQHHIFSFFAFKENRLPLFVKVRDTTQEPCGRLSFMKEPKSTRGLVHQAICNLNITLPIYAKESESDQEPEEEIGMTSEKNDETESTETSVLKSHLVNEVPVLASPDLLSEVSEMKQDLIKMTAILTTDVSDKAGSLKVKELAKAGEEEPGEPFEIVERVKEDLEKVNAILRSGTCMRDEGRARSSQSERELEEEWVIVSDEEIQEAKQHAPVEIDEHPCIEVRVDRETKAKVEKDSTGLVNYLTDDLNSYTSPHEKKPHTAPEKSGETSQASAVGKSSESNKGKATSAEEKQSAQKQLKPGLAIKKPVRRKLKEKQKQKEESSQSSEEKTELKKGSSEESVDEDRGLVPEPLPTAKATSPLIEETPIGSIKDKVKALQKRVEDEQKGRSKLPVRVKGKEDVPKRTTPRTHPAVSPSSKSSTSSKAERHSSLSSSAKPERHTPVSPSSKNEKLSPVSPSAKTERHSPVFSGKPEKHSPGSPSTKNERHSPVSSLKTERHTPGSPSGKTDKRPPVPSSGRTEKHPPVSPGKTEKHLPGSPSIRTPEKPAPGSATGKHEKHLPVSPGKTEKQPPISPTSKTERIEETMSVRELMKAFQSGQDPSKHKTGLFEHKSAKQKQPQDKSKSRVEKEKGHTVTQREVTQRETQRIESQTAKRGQRFQVSAATESRRFRSTTITVGLRMEDPVRERFERTPIIKTPEVVPSVAAEESHRGSEKIVDEQGDMDFQISPDRKTSTDFSEVIKQELEDNDKYQQFRLTEDTEKAQVHLDQVITSPFNTAFPLDYMKDEFLPALSLQSGALGGSSESLKQEVIAGSPCSSLMEGTPQISSEESYKHEGLAETPETSPESLSFSPKKSEEQIGEAKETTKVGTPTDIHSEKELPITNDITDSSQKQGAGVTRGSEPSTEHSQKEVTQDPHKDVCSKQDGCPESQSVSLASEVFTEKGSCGESQLPLVSSAFKTQSESETQESLTPSEVTKPFPPSDASVKTAEGTEPKPQGAIRSPQGLELPLPNRDSEVLSPMADESLAVSHKDSLEASPVLEDNSSHKTPDSLEPSPLKESPCRDSLESSPVEPKMKAGILPSHFPLPAAIAKTDLVAEVASMRSRLLRDPDGSAEDDSLEQTSLMESSGKSPLSPDTPSSEEVSYEVTPKPSDSSTPKPAVIHECAEEDDSENGEKKRFTPEEEMFKMVTKIKTFDELEQEAKQKRDYKKEPRQDGSSSASDPDADYSAEVNDEKQMAGTEGEGEVPVLVTSENRKVSSSSSESEPELTQLSKGADSGLLTEPVIRVQPPSPLPSSIDSNSSPEEATQFQPIVPKQYTFKMNEEIQEEPATSEDKDCKSHLAEDSQTHSADAADGSDGDLNRETTQPETCDGHGCETVSPSNSATPVSLGVQSPEHKDVDKPLAIDKDSLAHQDTCENDREEREFDPSGVESTQADLPNESSSLSSRCAIPEGNESAKEIASPSSPVKVEVTITDQALESMPEDCPIQDSSTTMQTERFAMDVPVSELAETDENSDPQIISPYENVPSSSFFSAEPSKIQTDTCHSTVVHSPEVYSVIIRSSPEDVVVTNSSNRTVSGEESHCESHDLETESEQKSALWAAQSDAPPLAVAPTASDAASVTGEQASKVIITKTDADADSWSEIREDDAAFEARVKEEEQKIFGLMVDRQSQGTTPDTTPARTPTEEGTPTSEQNPFLFQEGKLFEMTRSGAIDMTKRPYADESLHFFQIGQESNEEAISEDLKEGATGAEPPQTETTSESLELSEPKEAMDDEGELLPDDVSEEIEDLPASDANIDSQVIISASTETPTKEAVSTAVEEPPTTQRSDSLSTVKQTPRPAVPGPVGQLDFSPVTRSVYSGQDDESPESSPEEQKSVIEIPTAPVDNVPSAESKPQIPIRTLPTLVPAPPSAEDESAFSDDFPSSLDEDSKEGGAKPKSKIPVKAPTQRTEWQPSPTDIPLQKTAVPQGQETLSRAPDGRSKSESDASSLDAKTKCPVKARSYIETETESRERAEGFESESEDGATKPKLFASRLPVKSRSTSSSGRPGTSPTRESREHFFDLYRNSIEFFEEISDEASKLVDRLTQSEREQEPPSDDESSSALEVSVIESLPPVDIEHSAPEDIFDTRPIWDESIETMIERIPDENGHDRAEDPQDEQERMEERLAYIADHLGFSWTELARELDFTEEQIHQIRIENPNSLQDQSHALLKYWLERDGKHATDTILIECLTKINRMDIVHLLETNTEPLQERMGRSYAEIEQTITLDHSEGFSVLPDELCAAKEKKEQEASKESESSDHPPMVSEEDISVGYSTFQDCLPKTEGDSPAAALSPQMHQEPVQQDFSGKTQDQQEYYVTTPGAEVEDPQKATAVPDSLCKTPEDISTPPEGTKPCLQTPVTSERGSPIVQEPEEASEPKEESSPRKTSLVIVESTDDQSQVFERLDGDAAFQKGDDMPDIPPETVTEEEYVDENGHTVVKKVTRKIIRRYVSSDGTEKEEVTMQGMPQEPVNIEDGDNYSKVIKRVVLKSDTQQSEDNNE",
   "taxon_ids": [
    "10090"
   ],
   "accessions": [
    "Q8C8R3"
   ],
   "assertions": [
    {
     "assertion_id": "6c1efe2ee9c429e4832d89e49c0067c32db4c7400beed36961e6ba1035e4d66c:GO:0019901",
     "go_id": "GO:0019901",
     "go_name": "protein kinase binding",
     "go_definition": "\"Binding to a protein kinase, any enzyme that catalyzes the transfer of a phosphate group, usually from ATP, to a protein substrate.\" [GOC:jl]",
     "supporting_assertions": [
      {
       "accession": "Q8C8R3",
       "symbol": "Ank2",
       "relation": "enables",
       "go_id": "GO:0019901",
       "references": [
        "GO_REF:0000096"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "RGD:620156"
       ],
       "object_name": "Ankyrin-2",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000096.html"
       ]
      },
      {
       "accession": "Q8C8R3",
       "symbol": "Ank2",
       "relation": "enables",
       "go_id": "GO:0019901",
       "references": [
        "GO_REF:0000119"
       ],
       "evidence_code": "ISO",
       "with_from": [
        "UniProtKB:Q01484"
       ],
       "object_name": "Ankyrin-2",
       "taxon": [
        "taxon:10090"
       ],
       "annotation_date": "20260723",
       "assigned_by": "GO_Central",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000119.html",
        "https://www.uniprot.org/uniprotkb/Q01484/entry"
       ]
      }
     ]
    }
   ]
  },
  {
   "sequence_sha256": "928ccbd7cbbfc119d8ba7775d1093978032796836e86f47a5d0a392de6f28d5c",
   "sequence": "MGQIEWAMWANEQALASGLILITGGIVATAGRFTQWYFGAYSIVAGVFVCLLEYPRGKRKKGSTMERWGQKYMTAVVKLFGPFTRNYYVRAVLHLLLSVPAGFLLATILGTACLAIASGIYLLAAVRGEQWTPIEPKPRERPQIGGTIKQPPSNPPPRPPAEARKKPSEEEAAVAAGGPPGGPQVNPIPVTDEVV",
   "taxon_ids": [
    "9606"
   ],
   "accessions": [
    "P13498"
   ],
   "assertions": [
    {
     "assertion_id": "928ccbd7cbbfc119d8ba7775d1093978032796836e86f47a5d0a392de6f28d5c:GO:0020037",
     "go_id": "GO:0020037",
     "go_name": "heme binding",
     "go_definition": "\"Binding to a heme, a compound composed of iron complexed in a porphyrin (tetrapyrrole) ring.\" [GOC:ai]",
     "supporting_assertions": [
      {
       "accession": "P13498",
       "symbol": "CYBA",
       "relation": "enables",
       "go_id": "GO:0020037",
       "references": [
        "GO_REF:0000002"
       ],
       "evidence_code": "IEA",
       "with_from": [
        "InterPro:IPR007732"
       ],
       "object_name": "Cytochrome b-245 light chain",
       "taxon": [
        "taxon:9606"
       ],
       "annotation_date": "20260727",
       "assigned_by": "InterPro",
       "annotation_extension": "",
       "gene_product_form_id": "",
       "links": [
        "https://geneontology.org/GO_REF/0000002.html",
        "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR007732/"
       ]
      }
     ]
    }
   ]
  }
 ]
}
```
