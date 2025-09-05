comprehensive_cell_colors <- c(
  # Lineage colors
  "Epithelial" = "#0072B2",
  "Endothelial" = "#009E73",
  "Stromal" = "#D55E00",
  "Immune" = "#CC79A7",

  # Cell type colors
  "AT1" = "#AA0DFE",
  "AT2" = "#FC00BCFF",
  "Basal" = "#85660D",
  "Ciliated" = "#782AB6",
  "Club" = "#1C8356",
  "Goblet" = "#267358",
  "Adventitial.Fibroblasts" = "#F7E1A0",
  "Inflammatory.Fibroblasts" = "#E8C892",
  "Alveolar.Fibroblasts" = "#C4451C",
  "Peribronchial.Fibroblasts" = "#DFC27D",
  "Pericytes" = "#FF5000FF",
  "SMC" = "#F8A19F",
  "Aerocytes" = "#5EF1F2",
  "Arterial" = "#003380",
  "gCAP" = "#E0FF66",
  "Lymph" = "#8F7C00",
  "pulmonary-venous" = "#E896A9",
  "Venous-systemic" = "#4C005C",
  "Tip.Like.Cells" = "#45B7D1",
  "B.Cells" = "#90AD1C",
  "Mast.Cells" = "#B10DA1",
  "Neutrophils" = "#1CFFCE",
  "Plasma.Cells" = "#325A9B",
  "proliferating.cells" = "#8B008B",
  "Alveolar.Macrophages" = "#D16B16",
  "CD4.T.Cells" = "#C71585",
  "CD8.T.Cells" = "#FF1493",
  "CD8.T.Cytotoxic.Cells" = "#17BECF",
  "Classical.Monocytes" = "#9467BD",
  "Conventional.Dendritic.Cells.1" = "#8C564B",
  "Conventional.Dendritic.Cells.2" = "#E9967A",
  "Langerhans.Cells" = "#F781BF",
  "Megakaryocytes" = "#9932CC",
  "Monocyte.derived.Macrophages" = "#006400",
  "NK.Cells.Bright" = "#DAA520",
  "NK.Cells.Dim" = "#FFD700",
  "Non.Classical.Monocytes" = "#4B0082",
  "pDC" = "#00CED1",
  "Th17.Cells" = "#BDB76B",
  "Treg.Cells" = "#483D8B"
)
disease_colors <- c(
  "PAH" = "firebrick2",
  "CTRL" = "dodgerblue2"
)

comprehensive_cell_colors_v2 <- c(
  # Epithelial Lineage (Warmth, Earthiness, Richness)
  "AT2" = "#FF7F50",
  "AT1" = "#E65100",
  "Club" = "#D84315",
  "Basal" = "#A14A20",
  "Goblet" = "#8B0000",
  "Ciliated" = "#9B59B6",

  # Endothelial Lineage (Flow, Growth, Vibrant Distinctness)
  "pulmonary-venous" = "#2ECC71",
  "gCAP" = "#27AE60",
  "Arterial" = "#16A085",
  "Venous-systemic" = "#1ABC9C",
  "Lymph" = "#28B463",
  "Aerocytes" = "#6D9773",
  "Tip.Like.Cells" = "#33CCFF", # HIGHLIGHT

  # Immune Lineage (Dynamics, Coolness, Depth, Subtle Transitions)
  "CD4.T.Cells" = "#3498DB",
  "Th17.Cells" = "#2980B9",
  "CD8.T.Cells" = "#1F618D",
  "CD8.T.Cytotoxic.Cells" = "#154360",
  "Treg.Cells" = "#85C1E9",
  "B.Cells" = "#5DADE2",
  "Plasma.Cells" = "#AF7AC5",
  "NK.Cells.Dim" = "#7D3C98",
  "NK.Cells.Bright" = "#6C3483",
  "Conventional.Dendritic.Cells.1" = "#B0A9B0",
  "Conventional.Dendritic.Cells.2" = "#9D8AAD",
  "Langerhans.Cells" = "#76547C",
  "pDC" = "#8262B3",
  "Alveolar.Macrophages" = "#34495E",
  "Monocyte.derived.Macrophages" = "#4B677F",
  "Classical.Monocytes" = "#7F8C8D",
  "Non.Classical.Monocytes" = "#95A5A6",
  "Neutrophils" = "#D1D1D1",
  "Mast.Cells" = "#A9CCE3",
  "proliferating.cells" = "#728FCE",
  "Megakaryocytes" = "#A569BD",

  # Stromal Lineage (Foundation, Structure, Calm)
  "Pericytes" = "#7B3F00",
  "SMC" = "#8B4513",
  "Alveolar.Fibroblasts" = "#8B572A",
  "Peribronchial.Fibroblasts" = "#A0522D",
  "Adventitial.Fibroblasts" = "#B8860B",
  #"Adventitial.like.Fibroblasts" = "#CD7F32",
  "Inflammatory.Fibroblasts" = "#CD7F32"
  # Disease levels
  #"PAH" = "firebrick2",
  #"CTRL" = "dodgerblue2"

)

# subject colors for subsampled object

structure(c("rosybrown1", "azure", "darkolivegreen", "lightsalmon4",
            "lightblue2", "chocolate2", "orchid3", "royalblue4", "goldenrod3",
            "turquoise", "gold", "pink3", "thistle4", "lightcyan1", "tomato4",
            "darkseagreen", "seagreen1", "paleturquoise3", "aquamarine1",
            "lightsalmon3", "darkslateblue", "lightyellow1", "salmon4", "honeydew",
            "gold4", "brown", "navajowhite1", "steelblue3", "cyan", "darksalmon",
            "mediumorchid2", "lightgoldenrod4", "lemonchiffon", "deeppink4",
            "mediumpurple3", "royalblue1", "royalblue3"), names = c("CTRL_GLE822",
                                                                    "CTRL_GLE891", "CTRL_CTRL002", "CTRL_GLR635", "CTRL_CTRL007",
                                                                    "CTRL_GLE639", "CTRL_GLR609", "CTRL_GLR725", "PAH_EX13_2224",
                                                                    "PAH_EX14_7", "PAH_EX15_22900", "PAH_H20_16473", "PAH_H21_2610",
                                                                    "CTRL_CTR048", "CTRL_CTR068", "PAH_PAH004", "PAH_PAH005", "PAH_PAH006",
                                                                    "PAH_PAH007", "PAH_PAH008", "PAH_PAH009", "PAH_PAH010", "PAH_PAH011",
                                                                    "CTRL_CTR069", "CTRL_CTR070", "PAH_PAH012", "PAH_PAH013", "PAH_PAH014",
                                                                    "PAH_PAH015", "PAH_PAH016", "PAH_PAH017", "PAH_PAH018", "PAH_PAH019",
                                                                    "PAH_PAH020", "PAH_PAH021", "PAH_PAH022", "PAH_PAH023"))
