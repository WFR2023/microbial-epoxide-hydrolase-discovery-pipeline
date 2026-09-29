# Functional-context scoring system

The functional-context score is a heuristic prioritization aid applied to candidates that had already passed the preceding annotation- and sequence-based prioritization stages.

It is not a statistical probability, an experimentally calibrated confidence measure, or a direct estimate of enzymatic function.

## Component 1: prior candidate-priority class

- High-priority category: +3 points
- Moderate-priority category: +2 points
- Exploratory or other category: +1 point

## Component 2: within-species recurrence

The recurrence contribution is based on the number of analyzed assemblies of the corresponding species meeting the recurrence criterion:

- 3 or more assemblies: +2 points
- 1-2 assemblies: +1 point
- 0 assemblies: +0 points

Recurrence is restricted to the selected assembly set and must not be interpreted as population prevalence.

## Component 3: genomic-neighborhood categories

One point is added when at least one neighboring feature is assigned to each of the following prespecified categories:

- transport
- lipid metabolism
- oxidative stress
- detoxification/xenobiotic-associated functions
- virulence/host interaction
- regulation

Each category contributes a maximum of one point regardless of the number of neighboring features assigned to that category.

## Theoretical score range

The theoretical functional-context score ranges from 1 to 11 points.

The contextual-support classes used in the workflow are:

- score >= 8: high contextual support
- score 5-7: moderate contextual support
- score < 5: limited contextual support

These thresholds are heuristic and were not derived from statistical optimization, machine-learning training, or experimental calibration.

## Evidence kept separate from the numerical score

The following analytical layers are not independently added as numerical terms to the functional-context score:

- Pfam domain assignments
- motif observations
- phylogenetic placement
- structural-model confidence
- fpocket cavity characteristics
- molecular-docking results

BLASTp and annotation evidence contribute upstream to candidate prioritization but are not separately re-added as independent numerical terms in the functional-context score.

This separation avoids treating heterogeneous computational evidence layers as quantitatively equivalent.
