"""Classification of every candidate sentence: DISCARD (the quantifier is not a
claim about results) with reason, or CHECK (an absolute claim about results,
design or process that has a computation in s10_assolute_checks.py)."""
D_DEF = "definition of a method, measure or rule: the quantifier describes the procedure, not a result"
D_DES = "description of the simulation design or data: true by construction, not a result"
D_BG = "background or motivation, no claim about results"
D_COND = "the hit word is a condition ('without the check', 'random intercept only') or emphasis, not a universal quantifier; numbers in it are checked by the registry/number tracks"
D_SET = "the quantifier only names the subset (e.g. 'at most 25 clusters', 'all subsamples'); the numbers are point values checked by the number tracks"
D_HDR = "table header, row label or caption phrase; no claim"
D_MATH = "mathematical or software statement, true independently of the results"
D_PROC = "statement about process or availability (repository, versions, protocol), not checkable against the tables"
DISCARD = {
 "A0004":D_DEF,"A0005":D_BG,"A0006":D_BG,"A0009":D_DES,"A0010":D_DEF,"A0011":D_DES,"A0012":D_PROC,
 "A0013":D_DEF,"A0014":D_DEF,"A0017":D_DES,"A0019":D_DES,"A0023":D_DEF,"A0025":D_DEF,"A0026":D_DEF,
 "A0028":D_DEF,"A0029":D_MATH,"A0030":D_COND,"A0031":D_DEF,"A0033":D_MATH,"A0034":D_DEF,"A0035":D_DEF,
 "A0036":D_DEF,"A0040":D_DEF,"A0041":D_DEF,"A0042":D_DEF,"A0043":D_DEF,"A0044":D_DEF,
 "A0045":D_DEF,"A0046":D_DEF,"A0047":D_DEF,"A0048":D_DEF,"A0049":D_DEF,"A0050":D_DEF,"A0051":D_DEF,
 "A0053":D_DEF,"A0056":D_DEF,"A0057":D_DES,"A0058":D_DEF,"A0062":D_COND,"A0063":D_COND,"A0064":D_DES,
 "A0067":D_COND,"A0068":D_COND,"A0083":D_SET,"A0085":D_COND,"A0089":D_COND,"A0090":D_SET,
 "A0091":"'most failures' is not universal; the share is in Table S15 (number track)","A0092":D_COND,
 "A0097":"recommendation; 'needs no convergence check' is a property of GEE, not a result",
 "A0099":D_PROC,"A0100":D_COND,"A0101":D_DEF,"A0102":D_DEF,"A0107":D_MATH,"A0111":D_DES,"A0113":D_DES,
 "A0116":"repeats the GLMM-native truth claim checked under A0024","A0118":D_MATH,"A0119":D_DES,
 "A0120":D_COND,"A0121":D_BG,"A0129":D_PROC,"A0131":D_HDR,"A0132":D_DES,"A0134":D_DEF,
 "A0135":D_DES,"A0136":D_DES,"A0137":D_DES,"A0138":D_DES,"A0139":D_DES,"A0140":D_DES,"A0141":D_DES,
 "A0142":D_DES,"A0143":D_DES,"A0144":D_DEF,"A0146":D_HDR,"A0149":D_DEF,"A0151":D_DEF,"A0153":D_DEF,
 "A0154":D_HDR,"A0155":D_HDR,"A0156":D_HDR,"A0157":D_DEF,"A0160":D_DEF,"A0166":D_DEF,"A0167":D_DES,
 "A0168":D_DES,"A0169":D_DEF,"A0171":D_PROC,"A0172":D_DEF,"A0174":D_DEF,"A0175":D_DEF,"A0179":D_HDR,
 "A0183":D_HDR,"A0185":D_HDR,"A0187":D_HDR,"A0191":D_DEF,"A0192":D_COND,"A0195":D_HDR,"A0198":D_HDR,
 "A0200":D_HDR,"A0201":D_HDR,"A0202":D_HDR,"A0203":D_HDR,"A0204":D_HDR,"A0205":D_HDR,"A0206":D_DEF,
 "A0209":D_HDR,"A0211":D_PROC,"A0213":D_COND,"A0214":D_MATH,"A0216":D_COND,"A0220":D_HDR,"A0222":D_HDR,
 "A0224":D_HDR,"A0230":D_HDR,"A0233":D_HDR,"A0235":D_DEF,"A0237":D_DEF,"A0238":D_DEF,"A0239":D_DEF,
 "A0240":D_DEF,"A0241":D_DEF,"A0242":D_DEF,"A0245":D_SET,"A0246":D_COND,"A0247":D_COND,"A0248":D_SET,
 "A0249":D_DEF,"A0250":D_HDR,"A0252":D_DES,"A0253":D_DES,"A0254":D_COND,"A0257":D_COND,"A0259":D_HDR,
 "A0260":D_DES,"A0261":D_DEF,"A0262":D_DEF,"A0263":D_DEF,"A0264":D_COND,"A0265":D_HDR,"A0267":D_HDR,
 "A0268":D_DEF,"A0271":D_SET,"A0272":D_HDR,"A0274":"table row; its readings are checked under A0078/A0096",
 "A0275":D_DEF,"A0276":D_DEF,"A0278":D_DEF,"A0279":D_HDR,"A0280":D_MATH,"A0282":D_HDR,"A0284":D_HDR,
 "A0285":D_COND,"A0293":D_DEF,"A0296":D_HDR,"A0298":D_DES,"A0300":D_HDR,"A0304":D_HDR,"A0306":D_DEF,
 "A0307":D_HDR,"A0309":D_DEF,"A0310":D_DEF,"A0311":D_MATH,"A0312":D_COND,"A0313":D_COND,"A0314":D_COND,
 "A0317":D_MATH,"A0320":D_COND,"A0321":D_COND,"A0322":D_HDR,"A0323":D_HDR,"A0325":D_DEF,"A0327":D_DEF,
 "A0328":D_DEF,"A0329":D_HDR,"A0084":D_SET,"A0173":None,
}
DISCARD = {k:v for k,v in DISCARD.items() if v}
