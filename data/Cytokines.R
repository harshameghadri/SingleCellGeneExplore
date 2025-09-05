# First let's define comprehensive groups
cytokine_groups <- list(
  # Interleukins and receptors
  "Interleukins" = c(
    # Core interleukins
    "IL1A", "IL1B", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL8", "IL9", "IL10",
    "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL16", "IL17A", "IL18", "IL19",
    "IL20", "IL21", "IL22", "IL23A", "IL24", "IL25", "IL26", "IL27", "IL31", "IL32",
    "IL33", "IL34", "IL36A", "IL36B", "IL36G", "IL37", "IL38",
    # Receptors
    "IL1R1", "IL1R2", "IL2RA", "IL2RB", "IL2RG", "IL3RA", "IL4R", "IL5RA", "IL6R",
    "IL7R", "IL9R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1",
    "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL21R", "IL23R", "IL1RN"
  ),

  # Chemokines and receptors
  "Chemokines" = c(
    # CC chemokines
    "CCL1", "CCL2", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11",
    "CCL13", "CCL14", "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL24",
    # CXC chemokines
    "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCL10",
    "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17",
    # Receptors
    "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCR10",
    "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
    "ACKR1", "ACKR2", "ACKR3", "ACKR4"
  ),

  # TNF family and receptors
  "TNF_Family" = c(
    # Ligands
    "TNF", "LTA", "LTB", "TNFSF4", "CD40LG", "FASLG", "TNFSF10", "TNFSF11",
    "TNFSF12", "TNFSF13", "TNFSF14",
    # Receptors
    "TNFRSF1A", "TNFRSF1B", "CD40", "FAS", "TNFRSF4", "TNFRSF8", "TNFRSF9",
    "TNFRSF10A", "TNFRSF10B"
  ),

  # Growth factors and receptors
  "Growth_Factors" = c(
    # Factors
    "TGFB1", "TGFB2", "TGFB3", "PDGFA", "PDGFB", "VEGFA", "VEGFB", "VEGFC",
    "EGF", "FGF1", "FGF2", "HGF", "IGF1", "BMP2", "BMP4", "BMP7",
    # Receptors
    "TGFBR1", "TGFBR2", "PDGFRA", "PDGFRB", "FLT1", "KDR", "EGFR", "IGF1R"
  ),

  # Interferons and receptors
  "Interferons" = c(
    # Interferons
    "IFNA1", "IFNA2", "IFNB1", "IFNG", "IFNL1", "IFNL2",
    # Receptors
    "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2"
  ),

  # CSF family and receptors
  "CSF_Family" = c(
    "CSF1", "CSF2", "CSF3", "CSF1R", "CSF2RA", "CSF2RB", "CSF3R"
  )
)
