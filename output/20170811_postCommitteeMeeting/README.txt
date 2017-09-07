Slight reworking of the script structure for methylation analysis to make it more extensible for future analyses, incorporating other analysis methods and/or additional datasets. Contained here are outputs from this pass through the workflow.

New workflow: metadata.R --> process_idats.R --> preprocess_methylation.R --> substructure.R --> downstream EWAS, etc.

Multiple outputs exist for some scripts:
metadata -- added FRS calculation
substructure -- multiple re-runs in order to fix errors or incorporate new metadata into "nonMethData" object
testMRS_repCV -- try different train/test splits
