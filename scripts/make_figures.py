import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

v = pd.read_csv("results/variants/variants_table.tsv", sep="\t")
cov = pd.read_csv("results/alignment/depth_1kb.tsv", sep="\t",
                  header=None, names=["start", "depth"])
FLAGGED = {366125}  # strand-bias artifact (lacZ)

# Figure 1: variants by effect
counts = v["EFFECT"].value_counts()
fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(counts.index, counts.values, color="steelblue")
ax.invert_yaxis()
ax.set_xlabel("Number of variants")
ax.set_title("Filtered variants by predicted effect (n = %d)" % len(v))
fig.tight_layout()
fig.savefig("figures/variants_by_effect.png", dpi=200)

# Figure 2: coverage + variant map
colors = {"HIGH": "red", "MODERATE": "orange", "LOW": "gray", "MODIFIER": "lightgray"}
fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
a1.plot(cov["start"] / 1e6, cov["depth"], lw=0.6, color="steelblue")
a1.axvline(1.230767, ls="--", color="black", lw=0.8)
a1.text(1.25, 130, "umuD (deleted)", fontsize=8)
a1.set_ylim(0, 150)
a1.text(0.40, 120, "lac region spike (~24,500x, off scale)", fontsize=8)
a1.set_ylabel("Mean depth (1 kb windows)")
a1.set_title("Read depth and filtered variants along the genome")

for _, r in v.iterrows():
    marker = "x" if r["POS"] in FLAGGED else "o"
    a2.scatter(r["POS"] / 1e6, 0, c=colors.get(r["IMPACT"], "black"),
               marker=marker, s=30)
for _, r in v[v["GENE"].isin(["gyrA", "parC"])].iterrows():
    a2.annotate(r["GENE"], (r["POS"] / 1e6, 0), xytext=(r["POS"] / 1e6, 0.8),
                ha="center", fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.6))
from matplotlib.lines import Line2D
h = [Line2D([0],[0], marker="o", ls="", color=c, label=k) for k, c in colors.items()]
a2.legend(handles=h, loc="upper left", ncol=4, fontsize=7, frameon=False)
a2.set_ylim(-0.5, 1.6)
a2.set_yticks([])
a2.set_xlabel("Position (Mb)")
fig.tight_layout()
fig.savefig("figures/genome_overview.png", dpi=200)

# Figure 3: allele support (reads for the change)
v["alt_reads"] = v["AD"].str.split(",").str[1].astype(int)
v["ref_reads"] = v["AD"].str.split(",").str[0].astype(int)
v["alt_fraction"] = v["alt_reads"] / (v["alt_reads"] + v["ref_reads"])
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(v["alt_fraction"], bins=20, color="steelblue")
ax.set_xlabel("Fraction of reads supporting the variant")
ax.set_ylabel("Number of variants")
ax.set_title("Support for each variant")
fig.tight_layout()
fig.savefig("figures/variant_support.png", dpi=200)
print("done")
