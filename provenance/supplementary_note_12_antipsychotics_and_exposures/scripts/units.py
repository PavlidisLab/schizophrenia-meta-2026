"""The antipsychotic contrasts used (drug vs vehicle), grouped into units of animals / cultures that share controls.

One row per Gemma result set. `drugs` maps a regex over the factor value text to the drug label; every factor value
that matches enters the unit as its own contrast against the result set's baseline (vehicle). Units that share
samples (Kim 2018: the same design on array and RNA-seq) are one unit. Chosen 2026-09-26 from a Gemma 2 search for
haloperidol, clozapine, olanzapine, risperidone, aripiprazole, chlorpromazine, quetiapine, ziprasidone, sulpiride,
fluphenazine, lurasidone and 'antipsychotic' (data/gemma_search.json); excluded: chemogenetic clozapine-N-oxide / DREADD
designs, peripheral tissues, drug screens outside brain or neural cells, designs without a vehicle group.
"""
DRUGS = r"haloperidol|clozapine|olanzapine|risperidone|quetiapine|chlorpromazine|aripiprazole"

UNITS = [
    # unit, class, gemma id, GEO, species, region, result set, duration, dose/route (from the dataset description or factor)
    ("GSE6511_brain", "chronic", 382, "GSE6511", "mouse", "brain", 518925, "4 weeks", "per dataset description"),
    ("GSE6467_brain", "chronic", 383, "GSE6467", "mouse", "brain", 549533, "12 weeks", "per dataset description"),
    ("GSE45229_striatum", "chronic", 6233, "GSE45229", "mouse", "striatum", 477923, "chronic", "quetiapine 10, 100 mg/kg; haloperidol 0.3, 1 mg/kg"),
    ("GSE45229_frontal_cortex", "chronic", 6233, "GSE45229", "mouse", "frontal cortex", 477924, "chronic", "quetiapine 10, 100 mg/kg; haloperidol 0.3, 1 mg/kg"),
    ("GSE67723_striatum", "chronic", 13381, "GSE67723", "mouse", "striatum", 543220, "per description", "risperidone"),
    ("GSE67753_whole_brain", "chronic", 13699, "GSE67753", "mouse", "whole brain", 481003, "chronic", "haloperidol (Kim 2018, array)"),
    ("GSE67750_striatum", "chronic", 14610, "GSE67750", "mouse", "striatum", 530218, "chronic", "haloperidol (Kim 2018, RNA-seq; possibly the same animals as GSE67753, other tissue)"),
    ("GSE93918_PFC", "chronic", 14218, "GSE93918", "mouse", "prefrontal cortex", 521856, "chronic", "clozapine 10 mg/kg"),
    ("GSE119772_hypothalamus", "chronic", 14945, "GSE119772", "mouse", "hypothalamus", 519301, "per description", "olanzapine"),
    ("GSE180473_striatum", "chronic", 20169, "GSE180473", "mouse", "striatum", 506295, "28 days", "olanzapine 3.5 mg/kg/d; risperidone 1 mg/kg/d"),
    ("GSE268867_lateral_septum", "chronic", 37733, "GSE268867", "mouse", "lateral septum", 581698, "12 weeks", "olanzapine 50 mg/kg in chow (females)"),
    ("GSE2547_frontal_cortex", "chronic", 700, "GSE2547", "rat", "frontal cortex", 561416, "21 days", "olanzapine 2 mg/kg/d i.p."),
    ("GSE66275_3regions", "chronic", 13715, "GSE66275", "rat", "frontal cortex, hippocampus, striatum", 551546, "21 days", "haloperidol 0.25 mg/kg/d; risperidone 5 mg/kg/d"),
    ("GSE4031_cortex", "chronic", 33966, "GSE4031", "rat", "cerebral cortex", 575133, "per description", "haloperidol, clozapine at low / medium / high dose"),
    ("GSE227464_striatum", "chronic", 57920, "GSE227464", "rat", "striatum", 610253, "per description", "haloperidol"),
    ("GSE48954_striatum", "acute", 6714, "GSE48954", "mouse", "striatum", 553384, "1 to 8 h", "haloperidol, clozapine, risperidone (vehicle: tween 80)"),
    ("GSE59895_brain", "acute", 18717, "GSE59895", "rat", "brain", 571574, "0.25 to 5 d", "haloperidol 21 mg/kg (DrugMatrix)"),
    ("GSE163161_iN", "cells", 18011, "GSE163161", "human", "iPSC-derived excitatory neurons", 517518, "per description", "clozapine"),
    ("GSE36678_SKNSH", "cells", 9489, "GSE36678", "human", "SK-N-SH neuroblastoma", 543607, "6, 24, 48 h", "risperidone"),
    ("GSE197966_NT2N", "cells", 25070, "GSE197966", "human", "NT2-N neurons", 619816, "24 h", "quetiapine 50 uM"),
]
